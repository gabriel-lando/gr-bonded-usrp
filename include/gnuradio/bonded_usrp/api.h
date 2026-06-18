#pragma once

#include <gnuradio/attributes.h>

#ifdef gnuradio_bonded_usrp_EXPORTS
#define BONDED_USRP_API __GR_ATTR_EXPORT
#else
#define BONDED_USRP_API __GR_ATTR_IMPORT
#endif
