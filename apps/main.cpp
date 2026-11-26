#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "Vector.h"

using Vec3D = Math::Vec3D;


struct Volume{ Vec3D HalfExtent; };

struct SimulationState
{
    Volume m_Volume;
    double m_Temperatur;
    double BoltzmanFaktor;
    std::vector<Vec3D> m_EntityPositions; 
};

auto AcceptanceProbability( const SimulationState& State, const double EnergyDiffrence ) -> double
{
    double BoltzmanPropability = std::exp(-1 * State.BoltzmanFaktor * EnergyDiffrence);
    return std::min(1.0, BoltzmanPropability); 
}

auto isChangeAccepted( const SimulationState& State, double EnergyDiffrence ) -> bool 
{ 
    double prop = AcceptanceProbability(State, EnergyDiffrence);


    return true;
}

auto EnergyFunction( const SimulationState& State, const double epsilon, const double sigma ) -> double 
{
    constexpr double LennardJonesExponent = 6.0;
    const auto& Positions = State.m_EntityPositions;
    double ResultingEnergy = 0.0;
    for( size_t i = 0; i < Positions.size(); ++i )
    {
        for( size_t j = i+1; j < Positions.size(); ++j )
        {
            const double distance = (Positions[j] - Positions[i]).EukNorm();
            const double LJConstant = sigma / distance;
            ResultingEnergy += 4 * epsilon * ( std::pow(LJConstant, 2 * LennardJonesExponent) - std::pow(LJConstant, LennardJonesExponent) );
        }
    }

    return ResultingEnergy;
}

auto Simulate() -> double
{


    return 0.0;
}



int main( int argc, char* argv[] )
{

    std::cout << "Hello World" << "\n";

    return 0;
}
