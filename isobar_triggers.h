#ifndef ISOBAR_TRIGGERS_H
#define ISOBAR_TRIGGERS_H

#include <map>

class isobar_triggers {
    public:
    enum trigger_names {minbias, bht1_vpd30, bht1_vpd100};
    std::map<int, trigger_names> triggers {{600001, minbias},
                                           {600011, minbias},
                                           {600021, minbias},
                                           {600031, minbias},
                                           {600221, bht1_vpd30},
                                           {600231, bht1_vpd30},
                                           {600222, bht1_vpd100},
                                           {600232, bht1_vpd100}
                                          };
};



#endif // ISOBAR_TRIGGERS_H