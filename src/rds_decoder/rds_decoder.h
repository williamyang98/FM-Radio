#pragma once

#include <stdint.h>
#include <memory>
#include "utility/span.h"
#include "rds_constants.h"

class RDS_Decoder_Handler;
class LoggingBuffer;

class RDS_Decoder
{
private:
    RDS_Decoder_Handler* handler;
    std::unique_ptr<LoggingBuffer> logging_buffer; 
public:
    RDS_Decoder();
    ~RDS_Decoder();
    void SetHandler(RDS_Decoder_Handler* _handler) { handler = _handler; }
    void ProcessGroup(rds_group_t group);
private:
    bool OnGroupType(rds_group_t group, uint8_t code, bool version);
    bool OnGroup0A(rds_group_t group);
    bool OnGroup1A(rds_group_t group);
    bool OnGroup2A(rds_group_t group);
    bool OnGroup3A(rds_group_t group);
    bool OnGroup4A(rds_group_t group);
    bool OnGroup10A(rds_group_t group);
    bool OnGroup11A(rds_group_t group);
    bool OnGroup14A(rds_group_t group);
private:
    void PrintAltFreq(uint8_t x);
};