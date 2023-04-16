#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <filesystem>
#include <vector>
#include "utility/span.h"
#include "audio/frame.h"

namespace fs = std::filesystem;

class Audio_Scraper 
{
private:
    const fs::path dir;    
    FILE* fp_wav = NULL;
    int old_F_sample;
    int total_bytes_written = 0;

    std::vector<Frame<int16_t>> convert_buffer;
public:
    Audio_Scraper(const fs::path& _dir): dir(_dir) {
        fp_wav = NULL;
        total_bytes_written = 0;
    }
    ~Audio_Scraper();
    Audio_Scraper(const Audio_Scraper&) = delete;
    Audio_Scraper(Audio_Scraper&&) = delete;
    Audio_Scraper& operator=(const Audio_Scraper&) = delete;
    Audio_Scraper& operator=(Audio_Scraper&&) = delete;
    void on_audio_data(tcb::span<const Frame<float>> data, const int F_sample);
private:
    FILE* create_wav_file(const int F_sample);
    void update_wav_header(FILE* fp, const int nb_data_bytes);
    void close_wav_file(FILE* fp, const int nb_data_bytes);
};

class RDS_Scraper 
{
private:
    const fs::path dir;
    FILE* fp = NULL;
public:
    RDS_Scraper(const fs::path& _dir): dir(_dir) {}
    ~RDS_Scraper();
    void on_rds_bytes(tcb::span<const uint8_t> data);
};

class FM_Scraper 
{
private:
    const fs::path dir;
    Audio_Scraper audio_scraper;
    RDS_Scraper rds_scraper;
public:
    FM_Scraper(const char* _output_directory);
public:
    void on_audio_data(tcb::span<const Frame<float>> data, const int F_sample);
    void on_rds_bytes(tcb::span<const uint8_t> data);
};