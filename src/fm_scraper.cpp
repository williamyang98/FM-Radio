#include "./fm_scraper.h"
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <ctime>
#include <filesystem>
#include <limits.h>

#define LOG_MESSAGE(fmt, ...) fprintf(stderr, "[FM-Scraper]" fmt "\n", __VA_ARGS__)
#define LOG_ERROR(fmt, ...) fprintf(stderr, "error: [FM-Scraper]" fmt "\n", __VA_ARGS__)

const struct {
    bool is_stereo = true;
    int total_channels = 2;
    int bytes_per_sample = sizeof(int16_t);
} AUDIO_PARAMS;

template <typename ... U>
std::string create_string(U&& ... args) {
    constexpr size_t MAX_LENGTH = 512;
    char buffer[MAX_LENGTH];

    const int nb_written = snprintf(
        buffer, MAX_LENGTH-1, 
        std::forward<U>(args)...
    );

    if (nb_written <= 0) {
        return {};
    }

    std::string fmt_string { buffer, size_t(nb_written) };
    return fmt_string;
}

std::string get_current_time(void) {
    auto t = std::time(NULL);
    auto tm = *std::localtime(&t);

    return create_string(
        "%04d-%02d-%02dT%02d-%02d-%02d",
        tm.tm_year+1900, tm.tm_mon+1, tm.tm_mday,
        tm.tm_hour, tm.tm_min, tm.tm_sec
    );
}

Audio_Scraper::~Audio_Scraper() {
    if (fp_wav != NULL) {
        close_wav_file(fp_wav, total_bytes_written);
        fp_wav = NULL;
        total_bytes_written = 0;
    }
}

void Audio_Scraper::on_audio_data(tcb::span<const Frame<float>> data, const int F_sample) {
    if (old_F_sample != F_sample) {
        if (fp_wav != NULL) {
            close_wav_file(fp_wav, total_bytes_written);
            fp_wav = NULL;
            total_bytes_written = 0;
        }

        fp_wav = create_wav_file(F_sample);
        total_bytes_written = 0;
        old_F_sample = F_sample;
    }

    if (fp_wav == NULL) {
        return;
    }

    // grow the conversion buffer only
    const size_t N = data.size(); 
    if (convert_buffer.size() < N) {
        convert_buffer.resize(N);
    }

    // convert data
    constexpr float CONVERT_RESCALE = float(std::numeric_limits<int16_t>::max()) * 0.95f;
    for (size_t i = 0; i < N; i++) {
        convert_buffer[i] = Frame<int16_t>(data[i]*CONVERT_RESCALE);
    }

    const size_t nb_written = fwrite(convert_buffer.data(), sizeof(Frame<int16_t>), N, fp_wav);
    if (nb_written != N) {
        LOG_ERROR("[audio] Failed to write bytes %zu/%zu", nb_written, N);
    }
    total_bytes_written += int(nb_written);
    update_wav_header(fp_wav, total_bytes_written);
}

FILE* Audio_Scraper::create_wav_file(const int F_sample) {
    fs::create_directories(dir);
    auto time_string = get_current_time();
    auto filename = create_string("%s_audio.wav", time_string.c_str());
    auto filepath = dir / filename;
    auto filepath_str = filepath.string();

    FILE* fp = fopen(filepath_str.c_str(), "wb+");
    if (fp == NULL) {
        LOG_ERROR("[audio] Failed to open file %s", filepath_str.c_str());
        return fp;
    }

    LOG_MESSAGE("[audio] Opened file %s", filepath_str.c_str());

    // Source: http://soundfile.sapp.org/doc/WaveFormat/
    struct WavHeader {
        char     ChunkID[4];
        int32_t  ChunkSize;
        char     Format[4];
        // Subchunk 1 = format information
        char     Subchunk1ID[4];
        int32_t  Subchunk1Size;
        int16_t  AudioFormat;
        int16_t  NumChannels;
        int32_t  SampleRate;
        int32_t  ByteRate;
        int16_t  BlockAlign;
        int16_t  BitsPerSample;
        // Subchunk 2 = data 
        char     Subchunk2ID[4];
        int32_t  Subchunk2Size;
    } header;

    const int16_t NumChannels = AUDIO_PARAMS.total_channels;
    const int32_t BitsPerSample = AUDIO_PARAMS.bytes_per_sample * 8;
    const int32_t SampleRate = static_cast<int32_t>(F_sample);

    strncpy(header.ChunkID, "RIFF", 4);
    strncpy(header.Format, "WAVE", 4);
    strncpy(header.Subchunk1ID, "fmt ", 4);
    strncpy(header.Subchunk2ID, "data", 4);

    header.Subchunk1Size = 16;  // size of PCM format fields 
    header.AudioFormat = 1;     // Linear quantisation
    header.NumChannels = NumChannels;     
    header.SampleRate = SampleRate;
    header.BitsPerSample = BitsPerSample;
    header.ByteRate = header.SampleRate * header.NumChannels * header.BitsPerSample / 8;
    header.BlockAlign = header.NumChannels * header.BitsPerSample / 8;

    // We update these values when we close the file
    header.Subchunk2Size = 0;
    header.ChunkSize = 36 + header.Subchunk2Size; 

    fwrite(&header, sizeof(WavHeader), 1, fp);

    return fp;
}

void Audio_Scraper::update_wav_header(FILE* fp, const int nb_data_bytes) {
    const int32_t Subchunk2Size = nb_data_bytes;
    const int32_t ChunkSize = 36 + Subchunk2Size;

    // Source: http://soundfile.sapp.org/doc/WaveFormat/
    // Refer to offset of each field
    // ChunkSize
    fseek(fp, 4, SEEK_SET);
    fwrite(&ChunkSize, sizeof(int32_t), 1, fp);
    // Subchunk2Size
    fseek(fp, 40, SEEK_SET);
    fwrite(&Subchunk2Size, sizeof(int32_t), 1, fp);

    fseek(fp, 0, SEEK_END);
}

void Audio_Scraper::close_wav_file(FILE* fp, const int nb_data_bytes) {
    update_wav_header(fp, nb_data_bytes);
    fclose(fp);
}

RDS_Scraper::~RDS_Scraper() {
    if (fp != NULL) {
        fclose(fp);
    }
}

void RDS_Scraper::on_rds_bytes(tcb::span<const uint8_t> data) {
    if (fp == NULL) {
        fs::create_directories(dir);
        auto time_string = get_current_time();
        auto filename = create_string("%s_rds.bin", time_string.c_str());
        auto filepath = dir / filename;
        auto filepath_str = filepath.string();

        fp = fopen(filepath_str.c_str(), "wb+");
        if (fp == NULL) {
            LOG_ERROR("[rds] Failed to open file %s", filepath_str.c_str());
            return;
        }
        LOG_MESSAGE("[rds] Opened file %s", filepath_str.c_str());
    }

    const size_t N = data.size();
    const size_t nb_written = fwrite(data.data(), sizeof(uint8_t), N, fp);
    if (nb_written != N) {
        LOG_ERROR("[audio] Failed to write bytes %zu/%zu", nb_written, N);
    }
}

FM_Scraper::FM_Scraper(const char* _output_directory)
:   dir(_output_directory),
    audio_scraper(_output_directory),
    rds_scraper(_output_directory) 
{
    fs::create_directories(dir);
}

void FM_Scraper::on_audio_data(tcb::span<const Frame<float>> data, const int F_sample) {
    audio_scraper.on_audio_data(data, F_sample);
}

void FM_Scraper::on_rds_bytes(tcb::span<const uint8_t> data) {
    rds_scraper.on_rds_bytes(data);
}
