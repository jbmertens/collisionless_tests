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
  specs.initialization_type = SheetSimulation::initializationType::overdensity2d;
  SheetSimulation::Specs specs_alt = {0};

  specs_alt = specs;
  specs_alt.ns1 = 31; specs_alt.ns2 = 31;
  specs_alt.nx = 32; specs_alt.ny = 32;
  specs_alt.carriers_per_dx = 0;
  runSim(specs_alt, "sim_ns031_nx032_cpdx0", 40, false);

  specs_alt = specs;
  specs_alt.ns1 = 63; specs_alt.ns2 = 63;
  specs_alt.nx = 64; specs_alt.ny = 64;
  specs_alt.carriers_per_dx = 0;
  runSim(specs_alt, "sim_ns063_nx064_cpdx0", 40, false);

  specs_alt = specs;
  specs_alt.ns1 = 127; specs_alt.ns2 = 127;
  specs_alt.nx = 128; specs_alt.ny = 128;
  specs_alt.carriers_per_dx = 0;
  runSim(specs_alt, "sim_ns127_nx128_cpdx0", 40, false);
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

void runGaussianField()
{
  SheetSimulation::Specs specs = {0};
  specs.lx = 1.0; specs.ly = 1.0; specs.lz = 1.0;
  specs.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_ds;
  specs.deposit = SheetSimulation::depositScheme::PCS;
  specs.initialization_type = SheetSimulation::initializationType::gaussian_random;
  specs.dt = 0.02;

  // specs.nx = 8; specs.ny = 8; specs.nz = 8;
  // specs.ns1 = 4; specs.ns2 = 4; specs.ns3 = 4;
  // specs.carriers_per_dx = 2;
  // specs.carriers_per_dy = 2;
  // specs.carriers_per_dz = 2;
  // runSim(specs, "sim_GRF_nx08_ns04", 10, true);

  // specs.nx = 16; specs.ny = 16; specs.nz = 16;
  // specs.ns1 = 8; specs.ns2 = 8; specs.ns3 = 8;
  // specs.carriers_per_dx = 2;
  // specs.carriers_per_dy = 2;
  // specs.carriers_per_dz = 2;
  // runSim(specs, "sim_GRF_nx16_ns08", 10, true);

  // specs.nx = 32; specs.ny = 32; specs.nz = 32;
  // specs.ns1 = 16; specs.ns2 = 16; specs.ns3 = 16;
  // specs.carriers_per_dx = 2;
  // specs.carriers_per_dy = 2;
  // specs.carriers_per_dz = 2;
  // runSim(specs, "sim_GRF_nx32_ns16", 10, true);

  // specs.nx = 8; specs.ny = 8; specs.nz = 8;
  // specs.ns1 = 8; specs.ns2 = 8; specs.ns3 = 8;
  // specs.carriers_per_dx = 0;
  // specs.carriers_per_dy = 0;
  // specs.carriers_per_dz = 0;
  // runSim(specs, "sim_GRF_nx08_ns08", 500, true);

  // specs.nx = 16; specs.ny = 16; specs.nz = 16;
  // specs.ns1 = 16; specs.ns2 = 16; specs.ns3 = 16;
  // specs.carriers_per_dx = 0;
  // specs.carriers_per_dy = 0;
  // specs.carriers_per_dz = 0;
  // runSim(specs, "sim_GRF_nx16_ns16", 500, true);

  // specs.nx = 32; specs.ny = 32; specs.nz = 32;
  // specs.ns1 = 32; specs.ns2 = 32; specs.ns3 = 32;
  // specs.carriers_per_dx = 0;
  // specs.carriers_per_dy = 0;
  // specs.carriers_per_dz = 0;
  // runSim(specs, "sim_GRF_nx32_ns32", 500, true);

  // specs.nx = 64; specs.ny = 64; specs.nz = 64;
  // specs.ns1 = 64; specs.ns2 = 64; specs.ns3 = 64;
  // specs.carriers_per_dx = 0;
  // specs.carriers_per_dy = 0;
  // specs.carriers_per_dz = 0;
  // runSim(specs, "sim_GRF_nx64_ns64", 500, true);

  specs.nx = 8; specs.ny = 8; specs.nz = 8;
  specs.ns1 = 4; specs.ns2 = 4; specs.ns3 = 4;
  specs.carriers_per_dx = 2;
  specs.carriers_per_dy = 2;
  specs.carriers_per_dz = 2;
  runSim(specs, "sim_GRF_nx512_ns256_cpdx2", 500, true);
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
  specs.dt = 0.1;
  specs.deposit = SheetSimulation::depositScheme::PCS;
  specs.carrier_count_scheme = SheetSimulation::carrierCountScheme::per_dx;

  // runOverdensityTests(specs);
  // runUniformTests(specs);
  runGaussianField();

  return 0;
}
