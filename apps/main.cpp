#include <array>
#include <span>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <cfenv>

#include "Vector.h"

using Vec3D = Math::Vec3D;


constexpr double kb = 1.380649e-23;
constexpr size_t Dim = 3;

struct RandomNumbers
{
    std::mt19937 generator;
    std::uniform_real_distribution<> ProbabilityDistribution;
    std::uniform_int_distribution<> EntityDistribution;
    std::array<std::uniform_real_distribution<>, Dim> VolumeDistributions;
};

struct Volume{ Vec3D Extent; };

struct SimulationState
{
    Volume m_Volume;
    double m_Temperatur;
    double BoltzmanFaktor;
    std::vector<Vec3D> m_EntityPositions; 
};

class IPotential
{
public:
    virtual ~IPotential() = default;
    virtual auto getPotentialEnergy( const SimulationState& ) -> double = 0;
};

class LennardJonesPotential : public IPotential
{
    static constexpr double LJExponent = 6.0;
public:
    auto getPotentialEnergy( const SimulationState& State ) -> double override
    {
        const auto& Positions = State.m_EntityPositions;
        double ResultingEnergy = 0.0;
        for( size_t i = 0; i < Positions.size(); ++i )
        {
            for( size_t j = i+1; j < Positions.size(); ++j )
            {
                const double distance = (Positions[j] - Positions[i]).EukNorm();
                const double LJConstant = m_sigma / distance;
                ResultingEnergy += 4 * m_epsilon * ( std::pow(LJConstant, 2 * LJExponent) - std::pow(LJConstant, LJExponent) );
            }
        }

        return ResultingEnergy;
    }
public:
    LennardJonesPotential( double epsilon, double sigma ) : m_epsilon(epsilon), m_sigma(sigma) {}
private:
    double m_epsilon;
    double m_sigma;
};

auto AcceptanceProbability( const SimulationState& State, const double EnergyDiffrence ) -> double
{
    double BoltzmanPropability = std::exp(-1 * State.BoltzmanFaktor * EnergyDiffrence);
    return std::min(1.0, BoltzmanPropability); 
}

auto isChangeAccepted( const SimulationState& State, double EnergyDiffrence,  RandomNumbers& random ) -> bool 
{ 
    double prob = AcceptanceProbability(State, EnergyDiffrence);

    if( prob == 1.0 )
        return true;

    double rndm_Number = random.ProbabilityDistribution(random.generator);

    return rndm_Number <= prob;
}

void generateInitialKomposition( std::vector<Vec3D>& outPositions, RandomNumbers& random )
{
    for( auto& Entity : outPositions )
        for( size_t koord = 0; koord < Dim; ++koord )
            Entity[koord] = random.VolumeDistributions[koord](random.generator);
}

auto generateNewPosition( RandomNumbers& random ) -> Vec3D
{
    Vec3D res;
    for( size_t i = 0; i < Dim; ++i )
        res[i] = random.VolumeDistributions[i](random.generator);

    return res;
}

auto AveragePotentialEnergys( std::span<const double> Energys ) -> double
{
    double erg = 0.0;
    for( const auto& E : Energys )
        erg += E;

    return erg / (double) Energys.size();
}

//Fragen: Ab wann speichern wir die Energie werte. Wie Viele Energie werte speicher ich Welche epsilon und sigma werte ergebn eine schönes Ergebnis das man ausgeben kann
//welche und wie viele Startwerte geben wir in der Kommandozeile über
[[nodiscard]]
auto Simulate( Volume volume, double Temperatur, size_t NumEntitys, size_t seed, std::unique_ptr<IPotential> Potential ) -> double
{
    RandomNumbers random
        { 
            .generator=std::mt19937(seed), 
            .ProbabilityDistribution=std::uniform_real_distribution<>(0.0, 1.0),
            .EntityDistribution=std::uniform_int_distribution<>(0, static_cast<int>(NumEntitys)),
            .VolumeDistributions={
                std::uniform_real_distribution<> (0.0, volume.Extent[0]),
                std::uniform_real_distribution<> (0.0, volume.Extent[1]),
                std::uniform_real_distribution<> (0.0, volume.Extent[2])
            }
        };

    SimulationState State { .m_Volume=volume, .m_Temperatur=Temperatur, .BoltzmanFaktor=1/(Temperatur*kb), .m_EntityPositions=std::vector<Vec3D>(NumEntitys) };

    std::vector<double> ErgPotentialEnergy;

    generateInitialKomposition(State.m_EntityPositions, random);
    double KompositionEnergy = Potential->getPotentialEnergy(State);

    while(true)
    {
        auto randomEntity = static_cast<size_t>(random.EntityDistribution(random.generator));
        Vec3D PositionPrev = State.m_EntityPositions[randomEntity];

        State.m_EntityPositions[randomEntity] = generateNewPosition(random);
        double EnergyNew = Potential->getPotentialEnergy(State);

        double EnergyDelta = EnergyNew - KompositionEnergy;

        if( isChangeAccepted(State, EnergyDelta, random) )
        {
            KompositionEnergy = EnergyNew;
        }else {
            State.m_EntityPositions[randomEntity] = PositionPrev;
        }
    }

    return AveragePotentialEnergys(ErgPotentialEnergy);
}

int main( int argc, char* argv[] )
{
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    assert(argc >= 2);
    assert(std::atoi(argv[1]) > 0);

    auto NumEntitys = static_cast<size_t>(std::atoi(argv[1]));
    auto seed = static_cast<size_t>(std::atoi(argv[2]));

    if( seed == 0 ) 
    {
        std::random_device rd;
        seed = rd();
    }

    double Temperatur = 0.0;
    Volume volume { Vec3D{0.0, 0.0, 0.0} };
    std::unique_ptr<IPotential> Potential = std::make_unique<LennardJonesPotential>(0.0, 0.0);

    double erg = Simulate(volume, Temperatur, NumEntitys, seed, std::move(Potential));

    std::cout << "Avergae Potential Energy: "<< erg << "\n";

    return 0;
}
