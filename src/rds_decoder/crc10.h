#pragma once

#include <stdint.h>

// Calculate 10bit CRC from 26bit codeword
// 26bit -> 16bit data + 10bit crc
uint16_t CalculateCRC10(uint32_t x);

// Get error pattern as a 26bit bit string
// given the XOR difference between the received and transmitted CRC
// Return 0 if the syndrome doesn't correspond to a known error pattern
uint32_t GetCRCErrorFromSyndrome(uint16_t x);
