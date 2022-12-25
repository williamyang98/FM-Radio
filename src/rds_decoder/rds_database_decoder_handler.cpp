#include "rds_database_decoder_handler.h"
#include "rds_database.h"
#include <algorithm>
#include <stdio.h>

#define LOG_ERROR(...) fprintf(stderr, "ERROR: [rds_db_handler] " __VA_ARGS__)

RDS_Database_Decoder_Handler::RDS_Database_Decoder_Handler(RDS_Database& _db)
: db(_db)
{

}

void RDS_Database_Decoder_Handler::OnProgrammeIdentifier(uint16_t pi_code) {
    db.PI_code = pi_code;
}

void RDS_Database_Decoder_Handler::OnProgrammeType(uint8_t programme_type) {
    db.programme_type = programme_type;
}


void RDS_Database_Decoder_Handler::OnServiceName(char c, int index) {
    if (c == '\r') c = 0;
    db.service_name[index] = c;
}

void RDS_Database_Decoder_Handler::OnProgrammeTypeNameChange(uint8_t AB_flag) {
    if (AB_flag != AB_flag_programme_type_name) {
        std::memset(db.programme_type_name, 0, sizeof(db.programme_type_name));
    }
    AB_flag_programme_type_name = AB_flag;
}

void RDS_Database_Decoder_Handler::OnProgrammeTypeName(char c, int index) {
    if (c == '\r') c = 0;
    db.programme_type_name[index] = c;
}

void RDS_Database_Decoder_Handler::OnRadioTextChange(uint8_t AB_flag) {
    if (AB_flag != AB_flag_radio_text) {
        std::memset(db.radio_text, 0, sizeof(db.radio_text));
    }
    AB_flag_radio_text = AB_flag;
}

void RDS_Database_Decoder_Handler::OnRadioText(char c, int index) {
    if (c == '\r') c = 0;
    db.radio_text[index] = c;
}


void RDS_Database_Decoder_Handler::OnTrafficAnnouncement(bool traffic_announcement, bool traffic_programme) {
    // Clause 3.2.1.3: Traffic Programme (TP) and Traffic Announcement (TA) codes
    // Table 8: TP/TA
    uint8_t v = ((traffic_programme & 0b1) << 1) | (traffic_announcement & 0b1);
    switch (v) {
    case 0b00: 
        db.traffic_announcement = TrafficAnnouncement::NONE;
        break;
    case 0b01: 
        db.traffic_announcement = TrafficAnnouncement::EON_INFO;
        break;
    case 0b10: 
        db.traffic_announcement = TrafficAnnouncement::AWAIT_EON_ANNOUNCE;
        break;
    case 0b11: 
        db.traffic_announcement = TrafficAnnouncement::NOW_EON_ANNOUNCE;
        break;
    default:
        LOG_ERROR("Unknown TA/TP code: %u", v);
    }
}

void RDS_Database_Decoder_Handler::OnMusicSpeech(bool is_music) {
    db.is_music = is_music;
}


// Clause 3.2.1.5: Decoder Identification (DI) and Dynamic PTY Indicator (PTYI) codes
void RDS_Database_Decoder_Handler::OnDecoder_IsStereo(bool is_stereo) {
    db.is_stereo = is_stereo;
}

void RDS_Database_Decoder_Handler::OnDecoder_IsArtificalHead(bool is_artificial_head) {
    db.is_artificial_head = is_artificial_head;
}

void RDS_Database_Decoder_Handler::OnDecoder_IsCompressed(bool is_compressed) {
    db.is_compressed = is_compressed;
}

void RDS_Database_Decoder_Handler::OnDecoder_IsDynamicProgrammeType(bool is_dynamic) {
    db.is_dynamic_program_type = is_dynamic;
}


// Clause 3.2.1.6: Coding of Alternative Frequencies (AFs)
// Depending on the previous codes the meaning changes
void RDS_Database_Decoder_Handler::OnAlternativeFrequencyCode(uint8_t code, int index) {
    // TODO:
}


// Time and date
void RDS_Database_Decoder_Handler::OnDate(int day, int month, int year) {
    db.datetime.day = day;
    db.datetime.month = month;
    db.datetime.year = year;
}

void RDS_Database_Decoder_Handler::OnTime(uint8_t hour, uint8_t minute) {
    db.datetime.hour = hour;
    db.datetime.minute = minute;
}

void RDS_Database_Decoder_Handler::OnLocalTimeOffset(int8_t LTO) {
    db.local_time_offset = LTO;
}
