#include "rds_decoder.h"
#include "rds_constants.h"
#include "rds_programme_type_names.h"
#include "modified_julian_date.h"

#include "logging_buffer.h"
#include "rds_decoder_handler.h"
#include <stdio.h>

#if 1
#define APPEND_MESSAGE(...) logging_buffer->print(__VA_ARGS__)
#define LOG_MESSAGE(...) fprintf(stderr, "[rds_decoder] " __VA_ARGS__)
#else
#define APPEND_MESSAGE(...)
#define LOG_MESSAGE(...)
#endif

#define HANDLER if (handler) handler

// Clause 3.2.1.6.3: AF method A
void RDS_Decoder::PrintAltFreq(uint8_t x) {
    if (x == 0) {
        APPEND_MESSAGE("Unused");
        return;
    }

    // Table 11: Special meanings code table
    // Wait for this to occur to be able to contextualise following data
    if (x == 205) {
        APPEND_MESSAGE("Filler");
        return;
    }
    if ((x >= 224) && (x <= 249)) {
        const uint8_t total_AFs = x - 224;
        APPEND_MESSAGE("#AF%u", total_AFs);
        return;
    }  

    if (x == 250) {
        APPEND_MESSAGE("#LF/MF");
        return;
    }

    // TODO: VHF/LF/MF depends on the context bytes

    // Table 10: VHF code table (AF)
    if ((x >= 1) && (x <= 204)) {
        constexpr uint32_t base_freq = 87'500'000;
        constexpr uint32_t multiplier = 100'000;
        const uint32_t freq = base_freq + (uint32_t)x * multiplier;
        APPEND_MESSAGE("VHF=%.1fMHz", (float)freq * 1e-6f);
        return;
    }

    // Table 12: LF/MF code table - for ITU regions 1 and 3 (9 kHz spacing)
    if ((x >= 1) && (x <= 15)) {
        constexpr uint32_t base_freq = 153'000;
        constexpr uint32_t multiplier = 9'000;
        const uint32_t freq = base_freq + (uint32_t)(x-1)*multiplier;
        APPEND_MESSAGE("LF=%.1fkHz", (float)freq * 1e-3);
        return;
    }

    if ((x >= 16) && (x <= 135)) {
        constexpr uint32_t base_freq = 531'000;
        constexpr uint32_t multiplier = 9'000;
        const uint32_t freq = base_freq + (uint32_t)(x-16)*multiplier;
        APPEND_MESSAGE("MF=%.1fkHz", (float)freq * 1e-3);
        return;
    }

    APPEND_MESSAGE("Unassigned");
}

RDS_Decoder::RDS_Decoder() {
    handler = NULL;
    logging_buffer = std::make_unique<LoggingBuffer>();
}

RDS_Decoder::~RDS_Decoder() = default;

void RDS_Decoder::ProcessGroup(rds_group_t group) {
    // Begin text logging
    APPEND_MESSAGE("[group] [");
    for (int i = 0; i < RDS_PARAMS.nb_blocks_per_group; i++) {
        auto& block = group[i];
        if (block.is_valid) APPEND_MESSAGE("%04X", block.data);
        else                APPEND_MESSAGE("----");
        if (i != (RDS_PARAMS.nb_blocks_per_group-1)) APPEND_MESSAGE(" ");
        else                                         APPEND_MESSAGE("]");
    }

    // Clause 2.1: Baseband coding structure
    // Figure 9: Message format and addressing
    // Show the checksum status
    auto& block_A = group[0];
    auto& block_B = group[1];
    const uint16_t pi_code = block_A.data;
    const uint16_t descriptor = block_B.data;
    const uint8_t group_code    = (descriptor & 0b1111000000000000) >> 12; 
    const uint8_t version       = (descriptor & 0b0000100000000000) >> 11; 
    const uint8_t traffic_id    = (descriptor & 0b0000010000000000) >> 10;
    const uint8_t program_type  = (descriptor & 0b0000001111100000) >> 5;
    const uint8_t data          = (descriptor & 0b0000000000011111) >> 0;

    APPEND_MESSAGE(" ");
    if (block_A.is_valid) {
        HANDLER->OnProgrammeIdentifier(pi_code);
        APPEND_MESSAGE("PI=%04X, ", block_A.data);
    } else {
        APPEND_MESSAGE("         ");
    }

    if (block_B.is_valid) {
        APPEND_MESSAGE("Type %+2u%c, TP=%u, PTY=%+2u, ", 
            group_code, version ? 'B' : 'A',
            traffic_id,
            program_type);
        HANDLER->OnProgrammeType(program_type);
        OnGroupType(group, group_code, (bool)version);
    }

    auto buf = logging_buffer->c_str();
    LOG_MESSAGE("%.*s\n", (int)buf.size(), buf.data());
    logging_buffer->reset();
}

bool RDS_Decoder::OnGroupType(rds_group_t group, uint8_t code, bool version) {
    // Clause 3.1.3: Group types
    const bool is_version_A = !version;
    // Version A
    if (is_version_A) {
        switch (code) {
        case 0: return OnGroup0A(group);
        case 1: return OnGroup1A(group);
        case 2: return OnGroup2A(group);
        case 3: return OnGroup3A(group);
        case 4: return OnGroup4A(group);
        case 10: return OnGroup10A(group);
        case 11: return OnGroup11A(group);
        case 14: return OnGroup14A(group);
        default:
            APPEND_MESSAGE("Unsupported_Code");
            return false;
        }
    // Version B
    } else {
        switch (code) {
        default:
            APPEND_MESSAGE("Unsupported_Code");
            return false;
        }
    }

    // Table 3: Group types
    // Clause 3.2.1.3: Traffic Programme (TP) and Traffic Announcement (TA) codes
}

bool RDS_Decoder::OnGroup0A(rds_group_t group) {
    // Clause 3.1.5.1: Type 0 groups: Basic tuning and switching information
    // Figure 12: Basic tuning and switching information - Type 0A group
    auto& block_B = group[1];
    auto& block_C = group[2];
    auto& block_D = group[3];
    const bool has_block_C = block_C.is_valid && (block_C.block_type == BlockOffsetID::C);
    const bool has_block_D = block_D.is_valid && (block_D.block_type == BlockOffsetID::D);

    const uint8_t traffic_programme = (block_B.data & 0b0000'0100'0000'0000) >> 10;
    const uint8_t traffic_announce  = (block_B.data & 0b10000) >> 4;
    const uint8_t ms_flag           = (block_B.data & 0b01000) >> 3;
    const uint8_t decoder_bit       = (block_B.data & 0b00100) >> 2;
    const uint8_t segment_address   = (block_B.data & 0b00011) >> 0;

    // Clause 3.2.1.6.3: AF method A
    const uint8_t F0 = (block_C.data & 0xFF00) >> 8;
    const uint8_t F1 = (block_C.data & 0x00FF) >> 0;

    // Annex E: Character repertoires for Programme Service name, Programme Type Name, RadioText and alphanumeric Radio Paging
    constexpr int N_text = 2;
    char buf[N_text];
    buf[0] = has_block_D ? (char)((block_D.data & 0xFF00) >> 8) : '?';
    buf[1] = has_block_D ? (char)((block_D.data & 0x00FF) >> 0) : '?';

    HANDLER->OnMusicSpeech(ms_flag);
    HANDLER->OnTrafficAnnouncement(traffic_announce, traffic_programme);

    if (has_block_C) {
        const uint8_t index = 2*segment_address;
        HANDLER->OnAlternativeFrequencyCode(F0, index+0);
        HANDLER->OnAlternativeFrequencyCode(F1, index+1);
    }

    if (has_block_D) {
        const uint8_t index = 2*segment_address;
        HANDLER->OnServiceName(buf[0], index+0);
        HANDLER->OnServiceName(buf[1], index+1);
    }

    APPEND_MESSAGE("TA=%u, M/S=%u, decoder=%u, segment_address=%u, alt_freqs=[%03d,%03d] text='%.*s'",
        traffic_announce, 
        ms_flag, 
        decoder_bit, 
        segment_address,
        F0, F1,
        2, buf);

    // Clause 3.2.1.4: Music Speech (MS) switch code
    APPEND_MESSAGE(", ");
    APPEND_MESSAGE("M/S=%s", ms_flag ? "music" : "speech");

    // Clause 3.2.1.5: Decoder Identification (DI) and Dynamic PTY Indicator (PTYI) codes
    // Table 9: Bit d0 to d3 meanings
    APPEND_MESSAGE(", ");
    switch (segment_address) {
    case 0b00:  // d3
        HANDLER->OnDecoder_IsDynamicProgrammeType(decoder_bit);
        APPEND_MESSAGE("DI=%s", decoder_bit ? "dynamic_pty" : "static_pty");
        break;
    case 0b01:  // d2
        HANDLER->OnDecoder_IsCompressed(decoder_bit);
        APPEND_MESSAGE("DI=%s", decoder_bit ? "compressed" : "not_compressed");
        break;
    case 0b10:  // d1
        HANDLER->OnDecoder_IsArtificalHead(decoder_bit);
        APPEND_MESSAGE("DI=%s", decoder_bit ? "artificial_head" : "non_artificial_head");
        break;
    case 0b11:  // d0
        HANDLER->OnDecoder_IsStereo(decoder_bit);
        APPEND_MESSAGE("DI=%s", decoder_bit ? "stereo" : "mono");
        break;
    }

    APPEND_MESSAGE(", alt_freq=[");
    if (has_block_C) {
        PrintAltFreq(F0);
        APPEND_MESSAGE(",");
        PrintAltFreq(F1);
    } else {
        APPEND_MESSAGE("?,?");
    }
    APPEND_MESSAGE("]");

    return (has_block_C || has_block_D);
}

bool RDS_Decoder::OnGroup1A(rds_group_t group) {
    // Clause 3.1.5.2: Type 1 groups: Programme Item Number and slow labelling codes
    // Figure 14: Programme Item Number and slow labelling codes - Type 1A group
    auto& block_B = group[1];
    auto& block_C = group[2];
    auto& block_D = group[3];
    const bool has_block_C = block_C.is_valid && (block_C.block_type == BlockOffsetID::C);
    const bool has_block_D = block_D.is_valid && (block_D.block_type == BlockOffsetID::D);

    const uint8_t radio_paging_codes = (block_B.data & 0b0000'0000'0001'1111) >> 0;
    const uint8_t linkage_actuator   = (block_C.data & 0b1000'0000'0000'0000) >> 15;
    const uint8_t variant_code       = (block_C.data & 0b0111'0000'0000'0000) >> 12;
    const uint16_t data              = (block_C.data & 0b0000'1111'1111'1111) >> 0;
    const uint8_t day                = (block_D.data & 0b1111'1000'0000'0000) >> 11;
    const uint8_t hour               = (block_D.data & 0b0000'0111'1100'0000) >> 6;
    const uint8_t minute             = (block_D.data & 0b0000'0000'0011'1111) >> 0;

    APPEND_MESSAGE("radio_paging_code=%u, L/A=%u, variant=%u", 
        radio_paging_codes, linkage_actuator, variant_code);
    APPEND_MESSAGE(", ");

    switch (variant_code) {
    case 0b000:
        {
            const uint8_t paging                 = (data & 0b1111'0000'0000) >> 8;
            const uint16_t extended_country_code = (data & 0b0000'1111'1111) >> 0;
            APPEND_MESSAGE("paging=%u, ecc=%04X",
                paging, extended_country_code);
        }
        break;
    case 0b001: 
        APPEND_MESSAGE("tmc_id=%06X", data);   
        break;
    case 0b010: 
        APPEND_MESSAGE("paging_id=%06X", data);   
        break;
    case 0b011:
        APPEND_MESSAGE("language_code=%06X", data);   
        break;
    case 0b110:
        APPEND_MESSAGE("broadcast_use=%06X", data);   
        break;
    case 0b111:
        APPEND_MESSAGE("EWS_channel_id=%06X", data);   
        break;
    default:
        APPEND_MESSAGE("not_assigned_data=%06X", data);   
        break;
    }

    APPEND_MESSAGE(", ");
    APPEND_MESSAGE("day=%u, time=%02u:%02u", day, hour, minute);

    return (has_block_C || has_block_D);
}

bool RDS_Decoder::OnGroup2A(rds_group_t group) {
    // Clause 3.1.5.3: Type 2 groups: RadioText
    // Figure 16: RadioText - Type 2A group
    auto& block_B = group[1];
    auto& block_C = group[2];
    auto& block_D = group[3];
    const bool has_block_C = block_C.is_valid && (block_C.block_type == BlockOffsetID::C);
    const bool has_block_D = block_D.is_valid && (block_D.block_type == BlockOffsetID::D);

    const uint8_t AB_flag         = (block_B.data & 0b10000) >> 4;
    const uint8_t segment_address = (block_B.data & 0b01111) >> 0;

    constexpr int N_text = 4;
    char buf[N_text];
    buf[0] = has_block_C ? (char)((block_C.data & 0xFF00) >> 8) : '?';
    buf[1] = has_block_C ? (char)((block_C.data & 0x00FF) >> 0) : '?';
    buf[2] = has_block_D ? (char)((block_D.data & 0xFF00) >> 8) : '?';
    buf[3] = has_block_D ? (char)((block_D.data & 0x00FF) >> 0) : '?';

    const uint8_t index = segment_address*4;
    HANDLER->OnRadioTextChange(AB_flag);
    if (has_block_C) {
        HANDLER->OnRadioText(buf[0], index+0);
        HANDLER->OnRadioText(buf[1], index+1);
    }
    if (has_block_D) {
        HANDLER->OnRadioText(buf[2], index+2);
        HANDLER->OnRadioText(buf[3], index+3);
    }

    APPEND_MESSAGE("A/B=%u, segment_address=%+2u, text='%.*s'",
        AB_flag, segment_address,
        N_text, buf);

    return (has_block_C || has_block_D);
}

bool RDS_Decoder::OnGroup3A(rds_group_t group) {
    // Clause 3.1.5.4: Type 3A groups: Application identification for Open data
    // Figure 18: Application Identification for Open data - Type 3A group
    auto& block_B = group[1];
    auto& block_C = group[2];
    auto& block_D = group[3];
    const bool has_block_C = block_C.is_valid && (block_C.block_type == BlockOffsetID::C);
    const bool has_block_D = block_D.is_valid && (block_D.block_type == BlockOffsetID::D);

    const uint8_t app_code    = (block_B.data & 0b11111) >> 0;
    const uint8_t app_group   = (app_code & 0b11110) >> 1;
    const uint8_t app_version = (app_code & 0b00001) >> 0;

    const uint16_t message_bits = block_C.data;
    const uint16_t AID          = block_D.data;

    APPEND_MESSAGE("app_code=%u%c, message=%04X, AID=%04X", 
        app_group, app_version ? 'B' : 'A',
        message_bits,
        AID);

    return true;
}

bool RDS_Decoder::OnGroup4A(rds_group_t group) {
    // Clause 3.1.5.6: Type 4A groups Clock-time and date
    // Figure 20: Clock-time and date
    auto& block_B = group[1];
    auto& block_C = group[2];
    auto& block_D = group[3];
    const bool has_block_C = block_C.is_valid && (block_C.block_type == BlockOffsetID::C);
    const bool has_block_D = block_D.is_valid && (block_D.block_type == BlockOffsetID::D);

    const uint8_t rfu0     = ((block_B.data & 0b0000'0000'0001'1100) >> 2);
    const uint32_t MJD     = ((block_B.data & 0b0000'0000'0000'0011) << 15) |
                             ((block_C.data & 0b1111'1111'1111'1110) >> 1);
    const uint8_t hour     = ((block_C.data & 0b0000'0000'0000'0001) << 4) |
                             ((block_D.data & 0b1111'0000'0000'0000) >> 12);
    const uint8_t minute   = ((block_D.data & 0b0000'1111'1100'0000) >> 6);
    const uint8_t lto_sign = ((block_D.data & 0b0000'0000'0010'0000) >> 5);
    const uint8_t lto_val  = ((block_D.data & 0b0000'0000'0001'1111) >> 0);

    const int8_t LTO = (int8_t)lto_val * (lto_sign ? -1 : +1);

    int year, month, day;
    mjd_to_ymd((long)MJD, year, month, day);

    HANDLER->OnDate(day, month, year);
    HANDLER->OnTime(hour, minute);
    HANDLER->OnLocalTimeOffset(LTO);
    APPEND_MESSAGE("rfu0=%u, date=%02u/%02u/%04u, time=%02u:%02u, LTO=%d",
        rfu0,
        day, month, year,
        hour, minute,
        (int)LTO);

    return true;
}

bool RDS_Decoder::OnGroup10A(rds_group_t group) {
    // Clause 3.1.5.14: Type 10 groups: Programme Type Name (Group type 10A) and Open data (Group type 10B)
    // Figure 31: Programme Type Name PTYN - Type 10A group
    auto& block_B = group[1];
    auto& block_C = group[2];
    auto& block_D = group[3];
    const bool has_block_C = block_C.is_valid && (block_C.block_type == BlockOffsetID::C);
    const bool has_block_D = block_D.is_valid && (block_D.block_type == BlockOffsetID::D);

    const uint8_t AB_flag       = (block_B.data & 0b10000) >> 4;
    const uint8_t rfu0          = (block_B.data & 0b01110) >> 1;
    const uint8_t segment_addr  = (block_B.data & 0b00001) >> 0;

    constexpr int N_text = 4;
    char buf[N_text];
    buf[0] = has_block_C ? (char)((block_C.data & 0xFF00) >> 8) : '?';
    buf[1] = has_block_C ? (char)((block_C.data & 0x00FF) >> 0) : '?';
    buf[2] = has_block_D ? (char)((block_D.data & 0xFF00) >> 8) : '?';
    buf[3] = has_block_D ? (char)((block_D.data & 0x00FF) >> 0) : '?';

    const uint8_t index = 4*segment_addr;
    HANDLER->OnProgrammeTypeNameChange(AB_flag);
    if (has_block_C) {
        HANDLER->OnProgrammeTypeName(buf[0], index+0);
        HANDLER->OnProgrammeTypeName(buf[1], index+1);
    }
    if (has_block_D) {
        HANDLER->OnProgrammeTypeName(buf[2], index+2);
        HANDLER->OnProgrammeTypeName(buf[3], index+3);
    }

    APPEND_MESSAGE("A/B=%u, rfu0=%u, segment_addr=%u text='%.*s'",
        AB_flag, rfu0, segment_addr,
        N_text, buf);

    return true;
}

bool RDS_Decoder::OnGroup11A(rds_group_t group) {
    // Clause 3.1.5.15: Type 11 groups: Open Data Application
    // Figure 33: Open data - Type 11A and 11B groups
    // Clause 3.1.4 Open data channel / Applications Identification
    // TODO: This is not explicity specified in the standard
    APPEND_MESSAGE("TODO");
    return true;
}

bool RDS_Decoder::OnGroup14A(rds_group_t group) {
    // Clause 3.1.5.19: Type 14 groups: Enhanced Other Networks information
    // Figure 37: Enhanced Other Networks information - Type 14A groups
    auto& block_B = group[1];
    auto& block_C = group[2];
    auto& block_D = group[3];
    const bool has_block_C = block_C.is_valid && (block_C.block_type == BlockOffsetID::C);
    const bool has_block_D = block_D.is_valid && (block_D.block_type == BlockOffsetID::D);

    const uint8_t TP_on        = (block_B.data & 0b10000) >> 4;
    const uint8_t variant_code = (block_B.data & 0b01111) >> 0;
    const uint16_t data        =  block_C.data;
    const uint16_t PI_on       =  block_D.data;

    APPEND_MESSAGE("TP(on)=%u, variant=%u", TP_on, variant_code);
    APPEND_MESSAGE(", ");

    switch (variant_code) {
    case 0b0000:
    case 0b0001:
    case 0b0010:
    case 0b0011:
        {
            char buf[2];
            buf[0] = (char)((data & 0xFF00) >> 8);
            buf[1] = (char)((data & 0x00FF) >> 0);
            APPEND_MESSAGE("text='%.*s'", 2, buf);
        }
        break;
    case 0b0100:
        {
            const uint8_t F0 = (data & 0xFF00) >> 8;
            const uint8_t F1 = (data & 0x00FF) >> 0;
            APPEND_MESSAGE("AF(on)=[");
            PrintAltFreq(F0);
            APPEND_MESSAGE(",");
            PrintAltFreq(F1);
            APPEND_MESSAGE("]");
        }
        break;
    case 0b0101:
    case 0b0110:
    case 0b0111:
    case 0b1000:
        {
            // TODO: Unclear how the bitfield is structured from the figure
            // const uint8_t tuning_freq = data & 0b
            // const uint16_t mapped_fm_freq = data & 0b
            APPEND_MESSAGE("tuning_freq=?, mapped_fm_freq=?");
        }
        break;
    case 0b1001:
        {
            // TODO: Unclear how the bitfield is structured from the figure
            // const uint8_t tuning_freq = data & 0b
            // const uint16_t mapped_am_freq = data & 0b
            APPEND_MESSAGE("tuning_freq=?, mapped_am_freq=?");
        }
        break;
    case 0b1100:
        APPEND_MESSAGE("linkage_info=%04X", data);
        break;
    case 0b1101:
        {
            // TODO: Unclear how the bitfield is structured from the figure
            // const uint8_t PTY_on = data & 0b
            // const uint16_t rfu0 = data & 0b
            // const uint8_t TA = data & 0b
            APPEND_MESSAGE("bitfield_todo");
        }
        break;
    case 0b1110:
        APPEND_MESSAGE("PIN(on)=%04X", data);
        break;
    case 0b1111:
        APPEND_MESSAGE("reserved_broadcasters");
        break;
    default:
        APPEND_MESSAGE("Unallocated");
        break;
    }

    APPEND_MESSAGE(", ");
    APPEND_MESSAGE("PI(on)=%04X", PI_on);

    return true;
}
