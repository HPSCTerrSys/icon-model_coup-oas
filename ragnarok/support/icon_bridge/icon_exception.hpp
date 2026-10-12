// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_EXCEPTION_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_EXCEPTION_H_

#include <string>

namespace icon {

void message(const std::string& name, const std::string& text);
void finish(const std::string& name, const std::string& text);
void finish(const std::string& name);

}  // namespace icon

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_EXCEPTION_H_
