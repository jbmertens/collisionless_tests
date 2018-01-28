#define USE_OPENMP true

#include "sheet.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>
 
void runSim(SheetSimulation::Specs specs, std::string out_base_dir,
  int t_steps, bool writeslices)
{
  std::string mkdir_cmd ("mkdir -p " + out_base_dir + "; rm -f "
    + out_base_dir + "/*.gz");
  const int ret = system(mkdir_cmd.c_str());
  if(-1==ret) std::cout << "Unable to make dir: " + out_base_dir << std::endl;

  std::cout << std::setprecision(17);
  std::cout << "Initializing new sheet class with dir " << out_base_dir
    << "." << std::endl;
  SheetSimulation sheet (specs, SheetSimulation::Verbosity::none);
  sheet.writeInfo(out_base_dir);

  for(int t=0; t<=t_steps; ++t)
  {
    std::cout << "\rRunning step " << t << "..." << std::flush;
    sheet.writeStrips(out_base_dir + "/strips");
    sheet.writeConstraints(out_base_dir);
    if(writeslices && t % (t_steps/4) == 0)
    {
      std::string t_str = std::to_string(t);
      std::cout << std::endl << "Writing step " << t_str << std::endl;
      sheet.writeDensity(out_base_dir + "/step_" + t_str + "_density.dat");
      sheet.writePositions(out_base_dir + "/step_" + t_str + "_positions.dat");
    }
    sheet.RKStep();
  }
  std::cout << std::endl;
  sheet.printTimingInformation();
  std::cout << "Done!" << std::endl;
}


void runOverdensityTests(SheetSimulation::Specs specs)
{
  specs.initialization_type = SheetSimulation::initializationType::overdensity1d;
  SheetSimulation::Specs specs_alt = {0};

  specs_alt = specs;
  specs_alt.ns1 = 63;
  specs_alt.nx = 64;
  specs_alt.carriers_per_dx = 0;
  runSim(specs_alt, "sim_ns063_nx064_cpdx0", 400, false);

  specs_alt = specs;
  specs_alt.ns1 = 63;
  specs_alt.nx = 64;
  specs_alt.carriers_per_dx = 2;
  specs_alt.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_ds;
  runSim(specs_alt, "sim_ns063_nx064_cpds2", 400, false);


  specs_alt = specs;
  specs_alt.ns1 = 127;
  specs_alt.nx = 128;
  specs_alt.carriers_per_dx = 0;
  runSim(specs_alt, "sim_ns127_nx128_cpdx0", 400, false);

  specs_alt = specs;
  specs_alt.ns1 = 127;
  specs_alt.nx = 128;
  specs_alt.carriers_per_dx = 4;
  specs_alt.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_ds;
  runSim(specs_alt, "sim_ns127_nx128_cpds4", 400, false);


  specs_alt = specs;
  specs_alt.ns1 = 255;
  specs_alt.nx = 256;
  specs_alt.carriers_per_dx = 0;
  runSim(specs_alt, "sim_ns255_nx256_cpdx0", 400, false);

  specs_alt = specs;
  specs_alt.ns1 = 255;
  specs_alt.nx = 256;
  specs_alt.carriers_per_dx = 8;
  specs_alt.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_ds;
  runSim(specs_alt, "sim_ns255_nx256_cpds8", 400, false);
}


void runUniformTests(SheetSimulation::Specs specs)
{
  specs.initialization_type = SheetSimulation::initializationType::uniform1dv;
  SheetSimulation::Specs specs_alt = {0};

  specs_alt = specs;
  specs_alt.ns1 = 63;
  specs_alt.nx = 64;
  specs_alt.carriers_per_dx = 0;
  runSim(specs_alt, "sim_ns063_nx064_cpdx0", 400, false);

  specs_alt = specs;
  specs_alt.ns1 = 63;
  specs_alt.nx = 64;
  specs_alt.carriers_per_dx = 2;
  runSim(specs_alt, "sim_ns063_nx064_cpdx2", 400, false);

  specs_alt = specs;
  specs_alt.ns1 = 63;
  specs_alt.nx = 64;
  specs_alt.carriers_per_dx = 2;
  specs_alt.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_ds;
  runSim(specs_alt, "sim_ns063_nx064_cpds2", 400, false);


  specs_alt = specs;
  specs_alt.ns1 = 127;
  specs_alt.nx = 128;
  specs_alt.carriers_per_dx = 0;
  runSim(specs_alt, "sim_ns127_nx128_cpdx0", 400, false);

  specs_alt = specs;
  specs_alt.ns1 = 127;
  specs_alt.nx = 128;
  specs_alt.carriers_per_dx = 4;
  runSim(specs_alt, "sim_ns127_nx128_cpdx4", 400, false);

  specs_alt = specs;
  specs_alt.ns1 = 127;
  specs_alt.nx = 128;
  specs_alt.carriers_per_dx = 4;
  specs_alt.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_ds;
  runSim(specs_alt, "sim_ns127_nx128_cpds4", 400, false);


  specs_alt = specs;
  specs_alt.ns1 = 255;
  specs_alt.nx = 256;
  specs_alt.carriers_per_dx = 0;
  runSim(specs_alt, "sim_ns255_nx256_cpdx0", 400, false);

  specs_alt = specs;
  specs_alt.ns1 = 255;
  specs_alt.nx = 256;
  specs_alt.carriers_per_dx = 8;
  runSim(specs_alt, "sim_ns255_nx256_cpdx8", 400, false);

  specs_alt = specs;
  specs_alt.ns1 = 255;
  specs_alt.nx = 256;
  specs_alt.carriers_per_dx = 8;
  specs_alt.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_ds;
  runSim(specs_alt, "sim_ns255_nx256_cpds8", 400, false);

}

void runGaussianField(SheetSimulation::Specs specs)
{
  specs.nx = 32; specs.ny = 32; specs.nz = 32;
  specs.ns1 = 32; specs.ns2 = 32; specs.ns3 = 32;
  specs.carriers_per_dx = 2;
  specs.carriers_per_dy = 2;
  specs.carriers_per_dz = 2;
  specs.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_ds;
  specs.deposit = SheetSimulation::depositScheme::PCS;
  specs.initialization_type = SheetSimulation::initializationType::gaussian_random;

  SheetSimulation::Specs specs_alt = {0};

  specs_alt = specs;
  runSim(specs_alt, "sim_GRF", 100, true);
}


int main(int argc, char **argv)
{
  std::cout << "Beginning simulations." << std::endl;

  SheetSimulation::Specs specs = {0};
  specs.nx = 1; specs.ny = 1; specs.nz = 1;
  specs.ns1 = 1; specs.ns2 = 1; specs.ns3 = 1;
  specs.lx = 1.0; specs.ly = 1.0; specs.lz = 1.0;
  specs.carriers_per_dx = 4;
  specs.carriers_per_dy = 0;
  specs.carriers_per_dz = 0;
  specs.dt = 0.04;
  specs.deposit = SheetSimulation::depositScheme::PCS;
  specs.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_dx;

  // runOverdensityTests(specs);
  // runUniformTests(specs);
  runGaussianField(specs);

  return 0;
}
