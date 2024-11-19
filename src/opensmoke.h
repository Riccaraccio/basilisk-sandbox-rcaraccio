#include "OpenSMOKE_Interface.h"

#define OPENSMOKE 1

#pragma autolink -L$OPENSMOKE_INTERFACE/build -lopensmoke

char* kinfolder;

event cleanup (t = end)
{
  OpenSMOKE_Clean ();
}
