##### Prognostic variable

export DycoreVariables, Mass, Momentum, Energy
export Moisture, TotalMoisture, LiquidMoisture, IceMoisture
export Precipitation, Rain, Snow
export Tracers

abstract type DycoreVariables <: PrognosticVariable end
struct Mass <: DycoreVariables end
struct Momentum <: DycoreVariables end
struct Energy <: DycoreVariables end

abstract type Moisture <: PrognosticVariable end
struct TotalMoisture <: Moisture end
struct LiquidMoisture <: Moisture end
struct IceMoisture <: Moisture end

abstract type Precipitation <: PrognosticVariable end
struct Rain <: Precipitation end
struct Snow <: Precipitation end

struct Tracers <: PrognosticVariable end
