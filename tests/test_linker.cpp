#include "Teisko/BayerImage.hpp"
#include "Teisko/BayerInfo.hpp"
#include "Teisko/Chromaticity.hpp"
#include "Teisko/Color.hpp"
#include "Teisko/ColorCorrection.hpp"
#include "Teisko/LateralChromaticAberration.hpp"
#include "Teisko/LensShading.hpp"
#include "Teisko/MacbethDetector.hpp"
#include "Teisko/NoiseModel.hpp"
#include "Teisko/Preprocessing.hpp"
#include "Teisko/SpectralResponse.hpp"
#include "Teisko/Algorithm/Bit.hpp"
#include "Teisko/Algorithm/ConvexHull.hpp"
#include "Teisko/Algorithm/DelaunayTriangulation.hpp"
#include "Teisko/Algorithm/Functors.hpp"
#include "Teisko/Algorithm/Histogram.hpp"
#include "Teisko/Algorithm/Interpolate.hpp"
#include "Teisko/Algorithm/Iterators.hpp"
#include "Teisko/Algorithm/LinearSpace.hpp"
#include "Teisko/Algorithm/NelderMead.hpp"
#include "Teisko/Algorithm/PointXY.hpp"
#include "Teisko/Algorithm/Pow2.hpp"
#include "Teisko/Algorithm/ReduceTo.hpp"
#include "Teisko/Algorithm/StandardDeviation.hpp"
#include "Teisko/Algorithm/TrimmedMean.hpp"
#include "Teisko/Algorithm/VectorMedian.hpp"
#include "Teisko/Algorithm/VectorOperations.hpp"
#include "Teisko/Data/MunsellReflectance.hpp"
#include "Teisko/Image/Algorithms.hpp"
#include "Teisko/Image/API.hpp"
#include "Teisko/Image/Conversion.hpp"
#include "Teisko/Image/Point.hpp"
#include "Teisko/Image/Polyscale.hpp"
#include "Teisko/Image/RGB.hpp"
#include "Teisko/Image/Support.hpp"
#include "Teisko/Image/TIFF.hpp"

#include "catch.hpp"

TEST_CASE("Include everything")
{
}