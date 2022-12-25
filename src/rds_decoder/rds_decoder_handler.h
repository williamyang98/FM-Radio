#pragma once
#include <stdint.h>

class RDS_Decoder_Handler 
{
public:
    virtual void OnProgrammeIdentifier(uint16_t pi_code) = 0;
    virtual void OnProgrammeType(uint8_t programme_type) = 0;

    // Service names are expected to be static
    virtual void OnServiceName(char c, int index) = 0;
    virtual void OnProgrammeTypeName(char c, int index) = 0;
    virtual void OnRadioText(char c, int index) = 0;
    // Programme type name and radio text is expected to change
    // We do this by detecting a change in the AB_flag
    virtual void OnProgrammeTypeNameChange(uint8_t AB_flag) = 0;
    virtual void OnRadioTextChange(uint8_t AB_flag) = 0;

    virtual void OnTrafficAnnouncement(bool traffic_announcement, bool traffic_programme) = 0;
    virtual void OnMusicSpeech(bool is_music) = 0;

    // Clause 3.2.1.5: Decoder Identification (DI) and Dynamic PTY Indicator (PTYI) codes
    virtual void OnDecoder_IsStereo(bool is_stereo) = 0;
    virtual void OnDecoder_IsArtificalHead(bool is_artificial_head) = 0;
    virtual void OnDecoder_IsCompressed(bool is_compressed) = 0;
    virtual void OnDecoder_IsDynamicProgrammeType(bool is_dynamic) = 0;

    // Clause 3.2.1.6: Coding of Alternative Frequencies (AFs)
    // Depending on the previous codes the meaning changes
    virtual void OnAlternativeFrequencyCode(uint8_t code, int index) = 0;

    // Time and date
    virtual void OnDate(int day, int month, int year) = 0;
    virtual void OnTime(uint8_t hour, uint8_t minute) = 0;
    virtual void OnLocalTimeOffset(int8_t LTO) = 0;
};