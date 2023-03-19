#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <cstring>

typedef uint8_t programme_type_t; // 5bits
typedef uint32_t freq_t;          // Hz to MHz

struct AlternateFrequency {
    enum Type {
        VHF, MF, LF
    };
    Type type;
    freq_t freq;
};

enum TrafficAnnouncement {
    NONE,
    EON_INFO,
    AWAIT_EON_ANNOUNCE,
    NOW_EON_ANNOUNCE
};

struct RDS_Database 
{
    char service_name[8]{0};
    char programme_type_name[8]{0};
    char radio_text[64]{0};

    programme_type_t programme_type = 0;
    uint16_t PI_code = 0;

    bool is_stereo = false;
    bool is_music = false;
    bool is_artificial_head = false;
    bool is_compressed = false;
    bool is_dynamic_program_type = false;

    std::vector<AlternateFrequency> alt_freqs;

    struct datetime_t {
        int day = 0; 
        int month = 0; 
        int year = 0;
        uint8_t hour = 0;
        uint8_t minute = 0;
    } datetime;

    int8_t local_time_offset = 0;

    TrafficAnnouncement traffic_announcement = TrafficAnnouncement::NONE;

    void Reset() {
        memset(service_name, 0, sizeof(service_name));
        memset(programme_type_name, 0, sizeof(programme_type_name));
        memset(radio_text, 0, sizeof(radio_text));

        programme_type = 0;
        PI_code = 0;
        is_stereo = false;
        is_music = false;
        is_artificial_head = false;
        is_compressed = false;
        is_dynamic_program_type = false;

        datetime.day = 0;
        datetime.month = 0;
        datetime.year = 0;
        datetime.hour = 0;
        datetime.minute = 0;

        local_time_offset = 0;

        traffic_announcement = TrafficAnnouncement::NONE;

        alt_freqs.clear();
    }
};
