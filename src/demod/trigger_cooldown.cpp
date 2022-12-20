#include "trigger_cooldown.h"

// only propagate trigger signal if sample cooldown has passed
bool Trigger_Cooldown::on_trigger(bool trig) {
    if (trig && (N_remain == 0)) {
        N_remain = N_cooldown;
        return true;
    }
    if (N_remain > 0) {
        N_remain--;
    }
    return false;
}