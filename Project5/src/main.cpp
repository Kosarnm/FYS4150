#include "random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
	auto start = chrono::system_clock::now();
	int unitCellsPerDimention = 1;
	double initialTemperature = UnitConverter::temperatureFromSI(100.0); // measured in Kelvin
	double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms
	int numberOfTimeSteps = 10000;

	// If a first argument is provided, it is the number of unit cells
	if (numberOfArguments > 1) unitCellsPerDimention = atoi(argumentList[1]);
	// If a second argument is provided, it is the initial temperature (measured in kelvin)
	if (numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
	// If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
	if (numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));
	// If a third argument is provided, it is the number of time steps.
	if (numberOfArguments > 4) numberOfTimeSteps = atoi(argumentList[4]);


	double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.

	cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
	cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
	cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
	cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
	cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
	cout << "One unit of energy is " << UnitConverter::energyToSI(1.0) << " J" << endl;

	System system;
	system.createFCCLattice(unitCellsPerDimention, latticeConstant, initialTemperature);
	// system.create100Uniform(initialTemperature);

	system.potential().setEpsilon(1.0);
	system.potential().setSigma(3.405);

	//system.removeTotalMomentum();

	StatisticsSampler statisticsSampler;
	IO movie("movie.xyz"); // To write the state to file

	cout << setw(20) << "Timestep" <<
		setw(20) << "Time" <<
		setw(20) << "Temperature" <<
		setw(20) << "KineticEnergy" <<
		setw(20) << "PotentialEnergy" <<
		setw(20) << "TotalEnergy" << endl;
	for (int timestep = 0; timestep< numberOfTimeSteps; timestep++) {
		system.step(dt);
		statisticsSampler.sample(system);
		if (timestep % 100 == 0) {
			// Print the timestep every 100 timesteps
			cout << setw(20) << system.steps() <<
				setw(20) << system.time() <<
				setw(20) << statisticsSampler.temperature() <<
				setw(20) << statisticsSampler.kineticEnergy() <<
				setw(20) << statisticsSampler.potentialEnergy() <<
				setw(20) << statisticsSampler.totalEnergy() << endl;
			movie.saveState(system);
		}
	}

	movie.close();
	auto end = chrono::system_clock::now();
	chrono::duration<double> time_used = end - start;
	std::cout << time_used.count() << std::endl;

	return 0;
}