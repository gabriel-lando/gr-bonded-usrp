//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Export/visibility macro for the bonded engine compiled into the
// gr-bonded_usrp module. Replaces UHD_API now that the engine is no longer
// part of libuhd.
//

#pragma once

#if defined(_WIN32) || defined(__CYGWIN__)
#    define BONDED_API
#else
#    define BONDED_API __attribute__((visibility("default")))
#endif
