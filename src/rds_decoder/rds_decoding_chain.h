#pragma once

#include "rds_group_sync.h"
#include "rds_decoder.h"
#include "rds_database.h"
#include "rds_database_decoder_handler.h"

struct RDS_Decoding_Chain {
    RDS_Group_Sync group_sync;
    RDS_Decoder decoder;
    RDS_Database db;
    RDS_Database_Decoder_Handler db_handler;

    RDS_Decoding_Chain()
    : group_sync(), decoder(), db(),
      db_handler(db) 
    {
        group_sync.OnGroup().Attach([this](rds_group_t group) {
            decoder.ProcessGroup(group);
        });
        decoder.SetHandler(&db_handler);
    }

    void Process(tcb::span<const uint8_t> x) {
        group_sync.Process(x);
    }

    RDS_Decoding_Chain(RDS_Decoding_Chain&) = delete;
    RDS_Decoding_Chain(RDS_Decoding_Chain&&) = delete;
    RDS_Decoding_Chain& operator=(RDS_Decoding_Chain&) = delete;
    RDS_Decoding_Chain& operator=(RDS_Decoding_Chain&&) = delete;
};