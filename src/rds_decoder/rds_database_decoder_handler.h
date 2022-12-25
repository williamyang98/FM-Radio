#pragma once

#include "rds_decoder_handler.h"

struct RDS_Database;

class RDS_Database_Decoder_Handler: public RDS_Decoder_Handler
{
private:
    RDS_Database& db;
    uint8_t AB_flag_radio_text = 0b100;
    uint8_t AB_flag_programme_type_name = 0b100;
public:
    RDS_Database_Decoder_Handler(RDS_Database& _db);
    virtual void OnProgrammeIdentifier(uint16_t pi_code);
    virtual void OnProgrammeType(uint8_t programme_type);

    // Service names are expected to be static
    virtual void OnServiceName(char c, int index);
    virtual void OnProgrammeTypeName(char c, int index);
    virtual void OnRadioText(char c, int index);
    // Programme type name and radio text is expected to change
    // We do this by detecting a change in the AB_flag
    virtual void OnProgrammeTypeNameChange(uint8_t AB_flag);
    virtual void OnRadioTextChange(uint8_t AB_flag);

    virtual void OnTrafficAnnouncement(bool traffic_announcement, bool traffic_programme);
    virtual void OnMusicSpeech(bool is_music);

    // Clause 3.2.1.5: Decoder Identification (DI) and Dynamic PTY Indicator (PTYI) codes
    virtual void OnDecoder_IsStereo(bool is_stereo);
    virtual void OnDecoder_IsArtificalHead(bool is_artificial_head);
    virtual void OnDecoder_IsCompressed(bool is_compressed);
    virtual void OnDecoder_IsDynamicProgrammeType(bool is_dynamic);

    // Clause 3.2.1.6: Coding of Alternative Frequencies (AFs)
    // Depending on the previous codes the meaning changes
    virtual void OnAlternativeFrequencyCode(uint8_t code, int index);

    // Time and date
    virtual void OnDate(int day, int month, int year);
    virtual void OnTime(uint8_t hour, uint8_t minute);
    virtual void OnLocalTimeOffset(int8_t LTO);
};