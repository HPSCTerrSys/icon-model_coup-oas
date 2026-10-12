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
/// @brief read-in time-constant data fields without a time dimension
template <typename T>
void io::input_vector(netCDF::NcFile& datafile, std::vector<T>& v, const std::string& input, int ncells, int nlev) {
  netCDF::NcVar var;
  v.resize(ncells * nlev);
  //  access the input variable
  try {
    var = datafile.getVar(input);
  } catch (netCDF::exceptions::NcNotVar& e) {
    std::cerr << "FAILURE in accessing " << input << " (no time dimension): " << e.what() << std::endl;
    throw;
  }
  //  read-in input field values
  try {
    const std::vector<size_t> startp = {0, 0};
    const std::vector<size_t> count  = {static_cast<size_t>(nlev), static_cast<size_t>(ncells)};
    var.getVar(startp, count, v.data());
  } catch (netCDF::exceptions::NcNotVar& e) {
    std::cerr << "FAILURE in reading values from " << input << " (no time dimensions): " << e.what() << std::endl;
    throw;
  }
}

template <typename T>
void io::input_vector(netCDF::NcFile& datafile, std::vector<T>& v, const std::string& input, int ncells, int nlev,
                      int itime) {
  netCDF::NcVar att = datafile.getVar(input);
  try {
    v.resize(ncells * nlev);
    if (att.isNull()) {
      throw NC_ERR;
    }
    std::vector<size_t> startp = {static_cast<size_t>(itime), 0, 0};
    std::vector<size_t> count  = {1, static_cast<size_t>(nlev), static_cast<size_t>(ncells)};
    att.getVar(startp, count, v.data());
  } catch (netCDF::exceptions::NcException& e) {
    std::cerr << "FAILURE in reading " << input << ": " << e.what() << std::endl;
    throw;
  }
}

template <typename T, typename NCT>
void io::output_vector(netCDF::NcFile& datafile, std::vector<netCDF::NcDim>& dims, const std::string& output,
                       const std::vector<T>& v, int ncells, int nlev, int deflate_level) {
  // fortran:column major while c++ is row major
  NCT ncT;
  netCDF::NcVar var = datafile.addVar(output, ncT, dims);

  for (int i = 0; i < nlev; ++i) {
    var.putVar({static_cast<size_t>(i), 0}, {1, static_cast<size_t>(ncells)}, &v[i * ncells]);
  }
  if (deflate_level > 0) {
    var.setCompression(true, false, deflate_level);
  }
}

template <typename T, typename NCT>
void io::output_vector(netCDF::NcFile& datafile, std::vector<netCDF::NcDim>& dims, const std::string& output,
                       std::map<std::string, netCDF::NcVarAtt> varAttributes, const std::vector<T>& v, int ncells,
                       int nlev, int deflate_level) {
  // fortran:column major while c++ is row major
  NCT ncT;
  netCDF::NcVar var = datafile.addVar(output, ncT, dims);

  for (int i = 0; i < nlev; ++i) {
    var.putVar({static_cast<size_t>(i), 0}, {1, static_cast<size_t>(ncells)}, &v[i * ncells]);
  }
  if (deflate_level > 0) {
    var.setCompression(true, false, deflate_level);
  }
  // Add given attributes to the output variables (string, only)
  for (auto& attribute_name : {"standard_name", "long_name", "units", "coordinates", "CDI_grid_type"}) {
    auto attribute = varAttributes[attribute_name];

    // skip if attribute is not present
    if (attribute.isNull()) continue;

    std::string dataValues = "default";
    attribute.getValues(dataValues);

    var.putAtt(attribute.getName(), attribute.getType(), attribute.getAttLength(), dataValues.c_str());
  }
}

template <typename T>
void io::read_fields(const std::string& input_file, int itime, int& ncells, int& nlev, std::vector<T>& z,
                     std::vector<T>& t, std::vector<T>& p, std::vector<T>& rho, std::vector<T>& qv, std::vector<T>& qc,
                     std::vector<T>& qi, std::vector<T>& qr, std::vector<T>& qs, std::vector<T>& qg) {
  netCDF::NcFile datafile(input_file, netCDF::NcFile::read);

  //  read in the dimensions from the base variable: zg
  //  1st) vertical
  //  2nd) horizontal
  //
  auto baseDims = datafile.getVar(BASE_VAR).getDims();
  nlev          = baseDims[0].getSize();
  ncells        = baseDims[1].getSize();

  io::input_vector(datafile, z, "zg", ncells, nlev);
  io::input_vector(datafile, t, "ta", ncells, nlev, itime);

  io::input_vector(datafile, p, "pfull", ncells, nlev, itime);
  io::input_vector(datafile, rho, "rho", ncells, nlev, itime);
  io::input_vector(datafile, qv, "hus", ncells, nlev, itime);
  io::input_vector(datafile, qc, "clw", ncells, nlev, itime);
  io::input_vector(datafile, qi, "cli", ncells, nlev, itime);
  io::input_vector(datafile, qr, "qr", ncells, nlev, itime);
  io::input_vector(datafile, qs, "qs", ncells, nlev, itime);
  io::input_vector(datafile, qg, "qg", ncells, nlev, itime);

  datafile.close();
}

template <typename T, typename NCT>
void io::write_fields(const std::string& output_file, int ncells, int nlev, const std::vector<T>& t,
                      const std::vector<T>& qv, const std::vector<T>& qc, const std::vector<T>& qi,
                      const std::vector<T>& qr, const std::vector<T>& qs, const std::vector<T>& qg,
                      const std::vector<T>& prr_gsp, const std::vector<T>& pri_gsp, const std::vector<T>& prs_gsp,
                      const std::vector<T>& prg_gsp, const std::vector<T>& pflx, const std::vector<T>& pre_gsp) {
  netCDF::NcFile datafile(output_file, netCDF::NcFile::replace);
  netCDF::NcDim ncells_dim          = datafile.addDim("ncells", ncells);
  netCDF::NcDim nlev_dim            = datafile.addDim("height", nlev);
  std::vector<netCDF::NcDim> dims   = {nlev_dim, ncells_dim};
  int deflate_level                 = 0;
  int onelev                        = 1;
  netCDF::NcDim onelev_dim          = datafile.addDim("height1", onelev);
  std::vector<netCDF::NcDim> dims1d = {onelev_dim, ncells_dim};

  io::output_vector<T, NCT>(datafile, dims, "ta", t, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "hus", qv, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "clw", qc, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "cli", qi, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "qr", qr, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "qs", qs, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "qg", qg, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "pflx", pflx, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "prr_gsp", prr_gsp, ncells, onelev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "prs_gsp", prs_gsp, ncells, onelev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "pri_gsp", pri_gsp, ncells, onelev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "prg_gsp", prg_gsp, ncells, onelev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "pre_gsp", pre_gsp, ncells, onelev, deflate_level);

  datafile.close();
}

template <typename T>
static void copy_coordinate_variables_if_present(netCDF::NcFile& datafile, netCDF::NcFile& inputfile,
                                                 std::vector<std::string> coordinates) {
  for (auto& coordinate_name : coordinates) {
    auto coordinate = inputfile.getVar(coordinate_name);

    // skip if variable wasn't found
    if (coordinate.isNull()) continue;

    // copy possible new dimensions from input coordinates to the output
    // before adding the related data variables
    for (netCDF::NcDim& dim : coordinate.getDims()) {
      auto currentDims = datafile.getDims();
      // map.contains() would be better, but requires c++20
      if (auto search = currentDims.find(dim.getName()); search == currentDims.end())
        datafile.addDim(dim.getName(), dim.getSize());
    }
    auto var = datafile.addVar(coordinate.getName(), coordinate.getType(), coordinate.getDims());

    for (auto& [attribute_name, attribute] : coordinate.getAtts()) {
      std::string dataValues = "default";
      attribute.getValues(dataValues);
      var.putAtt(attribute.getName(), attribute.getType(), attribute.getAttLength(), dataValues.c_str());
    }

    // Calculate the total size of the variable
    auto dimensions = coordinate.getDims();
    int totalSize   = std::accumulate(dimensions.cbegin(), dimensions.cend(), 1,
                                      [](int acc, const netCDF::NcDim& dim) { return acc * dim.getSize(); });

    // Create a one-dimensional vector to store its values
    std::vector<T> oneDimensionalVariable(totalSize);

    // Copy the original values to the output
    coordinate.getVar(oneDimensionalVariable.data());
    var.putVar(oneDimensionalVariable.data());
  }
}

template <typename T, typename NCT>
void io::write_fields(const std::string& output_file, const std::string& input_file, int ncells, int nlev,
                      const std::vector<T>& t, const std::vector<T>& qv, const std::vector<T>& qc,
                      const std::vector<T>& qi, const std::vector<T>& qr, const std::vector<T>& qs,
                      const std::vector<T>& qg, const std::vector<T>& prr_gsp, const std::vector<T>& pri_gsp,
                      const std::vector<T>& prs_gsp, const std::vector<T>& prg_gsp, const std::vector<T>& pflx,
                      const std::vector<T>& pre_gsp) {
  netCDF::NcFile datafile(output_file, netCDF::NcFile::replace);
  netCDF::NcFile inputfile(input_file, netCDF::NcFile::read);
  auto baseDims            = inputfile.getVar(BASE_VAR).getDims();
  netCDF::NcDim nlev_dim   = datafile.addDim(baseDims[0].getName(), baseDims[0].getSize());
  netCDF::NcDim ncells_dim = datafile.addDim(baseDims[1].getName(), baseDims[1].getSize());

  copy_coordinate_variables_if_present<T>(
      datafile, inputfile, {"clon", "clon_bnds", "clat", "clat_bnds", baseDims[0].getName(), "height_bnds"});
  // height_bnds might have a different name

  std::vector<netCDF::NcDim> dims   = {nlev_dim, ncells_dim};
  int deflate_level                 = 0;
  int onelev                        = 1;
  netCDF::NcDim onelev_dim          = datafile.addDim("height1", onelev);
  std::vector<netCDF::NcDim> dims1d = {onelev_dim, ncells_dim};

  io::output_vector<T, NCT>(datafile, dims, "ta", inputfile.getVar("ta").getAtts(), t, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "hus", inputfile.getVar("hus").getAtts(), qv, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "clw", inputfile.getVar("clw").getAtts(), qc, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "cli", inputfile.getVar("cli").getAtts(), qi, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "qr", inputfile.getVar("qr").getAtts(), qr, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "qs", inputfile.getVar("qs").getAtts(), qs, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "qg", inputfile.getVar("qg").getAtts(), qg, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims, "pflx", pflx, ncells, nlev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "prr_gsp", prr_gsp, ncells, onelev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "prs_gsp", prs_gsp, ncells, onelev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "pri_gsp", pri_gsp, ncells, onelev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "prg_gsp", prg_gsp, ncells, onelev, deflate_level);
  io::output_vector<T, NCT>(datafile, dims1d, "pre_gsp", pre_gsp, ncells, onelev, deflate_level);

  inputfile.close();
  datafile.close();
}
