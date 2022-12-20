#pragma once

class Trigger_Cooldown
{
public:
    int N_cooldown = 1;
    int N_remain = 0;
    // only propagate trigger signal if sample cooldown has passed
    bool on_trigger(const bool trig);
};