//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
//    .--------------------------------------------------.
//    |  This file is part of NTCodeBase                 |
//    |  Created 2018 by NT (https://ttnghia.github.io)  |
//    '--------------------------------------------------'
//                            \o/
//                             |
//                            / |
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include <LibCommon/Utils/FileHelpers.h>
#include <LibCommon/Math/MathHelpers.h>
#include <LibCommon/ParallelHelpers/ParallelSTL.h>
#include <LibCommon/ParallelHelpers/Scheduler.h>
#include <LibParticle/ParticleHelpers.h>
#include <LibParticle/ParticleSerialization.h>

#include <Partio.h>

#include <fstream>
#include <cmath>
#include <vector>
#include <cassert>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase::ParticleHelpers {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
std::pair<VecX<N, Real_t>, VecX<N, Real_t>> getAABB(const StdVT_VecX<N, Real_t>& positions) {
    VecX<N, Real_t> bMin, bMax;
    ParallelSTL::min_max<N, Real_t>(positions, bMin, bMax);
    return std::make_pair(bMin, bMax);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> getCenter(const StdVT_VecX<N, Real_t>& positions) {
    auto [bMin, bMax] = getAABB(positions);
    return (bMin + bMax) * Real_t(0.5);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void compress(const StdVT_VecX<N, Real_t>& dvec, VecX<N, Real_t>& dMin, VecX<N, Real_t>& dMax, StdVT_UInt16& compressedData) {
    static_assert(std::is_floating_point_v<Real_t>);
    ////////////////////////////////////////////////////////////////////////////////
    ParallelSTL::min_max<N, Real_t>(dvec, dMin, dMax);
    auto diff = dMax - dMin;
    ////////////////////////////////////////////////////////////////////////////////
    // check if any dimension having the same value
    for(int j = 0; j < N; ++j) {
        if(std::abs(diff[j]) < Real_t(1e-20)) {
            diff[j] = Real_t(1); // to avoid divide by zero, setting to 1 is just fine
        }
    }
    ////////////////////////////////////////////////////////////////////////////////
    compressedData.resize(N * dvec.size());
    Scheduler::parallel_for(dvec.size(),
                            [&, diff](size_t i) {
                                const auto& vec = dvec[i];
                                for(int j = 0; j < N; ++j) {
                                    compressedData[i * N + j] = static_cast<UInt16>(std::numeric_limits<UInt16>::max() * ((vec[j] - dMin[j]) / diff[j]));
                                }
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void compress(const StdVT_VecX<N, Real_t>& dvec, DataBuffer& buffer, bool bWriteVectorSize /*= true*/) {
    static_assert(std::is_floating_point_v<Real_t>);
    ////////////////////////////////////////////////////////////////////////////////
    VecX<N, Real_t> dMin, dMax;
    StdVT_UInt16    compressedData;
    compress(dvec, dMin, dMax, compressedData);
    ////////////////////////////////////////////////////////////////////////////////
    // convert bmin and bmax to Vec3f
    VecX<N, float> dMinf, dMaxf;
    for(Int d = 0; d < N; ++d) {
        dMinf[d] = static_cast<float>(dMin[d]);
        dMaxf[d] = static_cast<float>(dMax[d]);
    }
    ////////////////////////////////////////////////////////////////////////////////
    buffer.clearBuffer();
    if(bWriteVectorSize) { buffer.append(static_cast<UInt>(dvec.size())); }
    buffer.append((const unsigned char*)glm::value_ptr(dMinf), sizeof(float) * N);
    buffer.append((const unsigned char*)glm::value_ptr(dMaxf), sizeof(float) * N);
    buffer.append((const unsigned char*)compressedData.data(), compressedData.size() * sizeof(UInt16));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void compress(const StdVT<MatXxX<N, Real_t>>& dvec, Real_t& dMin, Real_t& dMax, StdVT_UInt16& compressedData) {
    static_assert(std::is_floating_point_v<Real_t>);
    ////////////////////////////////////////////////////////////////////////////////
    Int NN = N * N;
    ParallelSTL::min_max<N, Real_t>(dvec, dMin, dMax);
    Real_t diff = dMax - dMin;
    ////////////////////////////////////////////////////////////////////////////////
    // check if any dimension having the same value
    if(std::abs(diff) < Real_t(1e-20)) {
        diff = Real_t(1); // to avoid divide by zero, setting to 1 is just fine
    }
    ////////////////////////////////////////////////////////////////////////////////
    compressedData.resize(NN * dvec.size());
    Scheduler::parallel_for(dvec.size(),
                            [&](size_t i) {
                                MatXxX<N, Real_t> mat = dvec[i];
                                const Real_t* mdata   = glm::value_ptr(mat);
                                for(int j = 0; j < NN; ++j) {
                                    compressedData[i * NN + j] = static_cast<UInt16>(std::numeric_limits<UInt16>::max() * ((mdata[j] - dMin) / diff));
                                }
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void compress(const StdVT<MatXxX<N, Real_t>>& dvec, DataBuffer& buffer, bool bWriteVectorSize /*= true*/) {
    static_assert(std::is_floating_point_v<Real_t>);
    ////////////////////////////////////////////////////////////////////////////////
    Real_t       dMin, dMax;
    StdVT_UInt16 compressedData;
    compress(dvec, dMin, dMax, compressedData);

    float dMinf = static_cast<float>(dMin);
    float dMaxf = static_cast<float>(dMax);
    buffer.clearBuffer();
    if(bWriteVectorSize) { buffer.append(static_cast<UInt>(dvec.size())); }
    buffer.append(&dMinf,                                      sizeof(float));
    buffer.append(&dMaxf,                                      sizeof(float));
    buffer.append((const unsigned char*)compressedData.data(), compressedData.size() * sizeof(UInt16));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void compress(const StdVT<Real_t>& dvec, Real_t& dMin, Real_t& dMax, StdVT_UInt16& compressedData) {
    static_assert(std::is_floating_point_v<Real_t>);
    ////////////////////////////////////////////////////////////////////////////////
    ParallelSTL::min_max<Real_t>(dvec, dMin, dMax);
    Real_t diff = dMax - dMin;
    ////////////////////////////////////////////////////////////////////////////////
    // check if any dimension having the same value
    if(std::abs(diff) < Real_t(1e-20)) {
        diff = Real_t(1); // to avoid divide by zero, setting to 1 is just fine
    }
    ////////////////////////////////////////////////////////////////////////////////
    compressedData.resize(dvec.size());
    Scheduler::parallel_for(dvec.size(),
                            [&](size_t i) {
                                compressedData[i] = static_cast<UInt16>(std::numeric_limits<UInt16>::max() * ((dvec[i] - dMin) / diff));
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void compress(const StdVT<Real_t>& dvec, DataBuffer& buffer, bool bWriteVectorSize /*= true*/) {
    static_assert(std::is_floating_point_v<Real_t>);
    ////////////////////////////////////////////////////////////////////////////////
    Real_t       dMin, dMax;
    StdVT_UInt16 compressedData;
    compress(dvec, dMin, dMax, compressedData);

    // convert bmin and bmax to Vec3f
    float dMinf = static_cast<float>(dMin);
    float dMaxf = static_cast<float>(dMax);
    buffer.clearBuffer();
    if(bWriteVectorSize) { buffer.append(static_cast<UInt>(dvec.size())); }
    buffer.append((const unsigned char*)&dMinf,                sizeof(float));
    buffer.append((const unsigned char*)&dMaxf,                sizeof(float));
    buffer.append((const unsigned char*)compressedData.data(), compressedData.size() * sizeof(UInt16));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void compress(const StdVT<StdVT<Real_t>>& dvec, StdVT<Real_t>& dMin, StdVT<Real_t>& dMax, StdVT<StdVT_UInt16>& compressedData) {
    static_assert(std::is_floating_point_v<Real_t>);
    __NT_REQUIRE(dvec.size() == dMin.size() && dvec.size() == dMax.size());
    ////////////////////////////////////////////////////////////////////////////////
    compressedData.resize(dvec.size());
    Scheduler::parallel_for(dvec.size(), [&](size_t i) { compress(dvec[i], dMin[i], dMax[i], compressedData[i]); });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void compress(const StdVT<StdVT<Real_t>>& dvec, DataBuffer& buffer, bool bWriteVectorSize /*= true*/) {
    static_assert(std::is_floating_point_v<Real_t>);
    StdVT<Real_t>       dMin, dMax;
    StdVT<StdVT_UInt16> compressedData;
    compress(dvec, dMin, dMax, compressedData);

    StdVT<float> dMinf(dvec.size());
    StdVT<float> dMaxf(dvec.size());

    for(size_t i = 0; i < dvec.size(); ++i) {
        dMinf[i] = static_cast<float>(dMin[i]);
        dMaxf[i] = static_cast<float>(dMax[i]);
    }

    buffer.clearBuffer();
    if(bWriteVectorSize) { buffer.append(static_cast<UInt>(dvec.size())); }
    for(size_t i = 0; i < dvec.size(); ++i) {
        buffer.append(static_cast<UInt>(compressedData[i].size()));
        buffer.append((const unsigned char*)&dMinf[i],                sizeof(float));
        buffer.append((const unsigned char*)&dMaxf[i],                sizeof(float));
        buffer.append((const unsigned char*)compressedData[i].data(), compressedData[i].size() * sizeof(UInt16));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void decompress(StdVT_VecX<N, Real_t>& dvec, const VecX<N, Real_t>& dMin, const VecX<N, Real_t>& dMax, const StdVT_UInt16& compressedData) {
    static_assert(std::is_floating_point_v<Real_t>);
    const VecX<N, Real_t> diff = dMax - dMin;
    __NT_REQUIRE((compressedData.size() / N) * N == compressedData.size());
    dvec.resize(compressedData.size() / N);
    Scheduler::parallel_for(dvec.size(),
                            [&](size_t i) {
                                VecX<N, Real_t> vec;
                                for(int j = 0; j < N; ++j) {
                                    vec[j] = static_cast<typename VecX<N, Real_t>::value_type>(compressedData[i * N + j]) * diff[j] /
                                             static_cast<typename VecX<N, Real_t>::value_type>(std::numeric_limits<UInt16>::max()) + dMin[j];
                                }
                                dvec[i] = vec;
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void decompress(StdVT_VecX<N, Real_t>& dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/) {
    static_assert(std::is_floating_point_v<Real_t>);
    VecX<N, float> dMinf, dMaxf;

    UInt64 segmentStart = 0;
    UInt64 segmentSize;

    if(nParticles == 0) {
        segmentSize = sizeof(UInt);
        memcpy(&nParticles, &buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;
    }

    segmentSize = sizeof(float) * N;
    memcpy(glm::value_ptr(dMinf), &buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;
    memcpy(glm::value_ptr(dMaxf), &buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;

    segmentSize = nParticles * N * sizeof(UInt16);
    __NT_REQUIRE(segmentStart + segmentSize == buffer.size());
    StdVT_UInt16 compressedData(nParticles * N);
    memcpy(compressedData.data(), &buffer.data()[segmentStart], segmentSize);

    VecX<N, Real_t> dMin, dMax;
    for(Int d = 0; d < N; ++d) {
        dMin[d] = static_cast<Real_t>(dMinf[d]);
        dMax[d] = static_cast<Real_t>(dMaxf[d]);
    }
    decompress(dvec, dMin, dMax, compressedData);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void decompress(StdVT<MatXxX<N, Real_t>>& dvec, Real_t dMin, Real_t dMax, const StdVT_UInt16& compressedData) {
    static_assert(std::is_floating_point_v<Real_t>);
    Int          NN   = N * N;
    const Real_t diff = dMax - dMin;
    __NT_REQUIRE((compressedData.size() / NN) * NN == compressedData.size());

    dvec.resize(compressedData.size() / NN);
    Scheduler::parallel_for(dvec.size(),
                            [&](size_t i) {
                                MatXxX<N, Real_t> mat;
                                Real_t* mdata = glm::value_ptr(mat);

                                for(int j = 0; j < NN; ++j) {
                                    mdata[j] = static_cast<typename VecX<N, Real_t>::value_type>(compressedData[i * NN + j]) * diff /
                                               static_cast<typename VecX<N, Real_t>::value_type>(std::numeric_limits<UInt16>::max()) + dMin;
                                }
                                dvec[i] = mat;
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void decompress(StdVT<MatXxX<N, Real_t>>& dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/) {
    static_assert(std::is_floating_point_v<Real_t>);
    float  dMinf, dMaxf;
    UInt64 segmentStart = 0;
    UInt64 segmentSize;

    if(nParticles == 0) {
        segmentSize = sizeof(UInt);
        memcpy(&nParticles, &buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;
    }

    segmentSize = sizeof(float);
    memcpy(&dMinf, &buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;
    memcpy(&dMaxf, &buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;

    segmentSize = nParticles * N * N * sizeof(UInt16);
    __NT_REQUIRE(segmentStart + segmentSize == buffer.size());
    StdVT_UInt16 compressedData(nParticles * N * N);
    memcpy(compressedData.data(), &buffer.data()[segmentStart], segmentSize);

    Real_t dMin = static_cast<Real_t>(dMinf);
    Real_t dMax = static_cast<Real_t>(dMaxf);
    decompress(dvec, dMin, dMax, compressedData);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void decompress(StdVT<Real_t>& dvec, Real_t dMin, Real_t dMax, const StdVT_UInt16& compressedData) {
    static_assert(std::is_floating_point_v<Real_t>);
    const Real_t diff = dMax - dMin;
    dvec.resize(compressedData.size());
    Scheduler::parallel_for(dvec.size(),
                            [&](size_t i) {
                                dvec[i] = static_cast<Real_t>(compressedData[i]) * diff / static_cast<Real_t>(std::numeric_limits<UInt16>::max()) + dMin;
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void decompress(StdVT<Real_t>& dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/) {
    static_assert(std::is_floating_point_v<Real_t>);
    float  dMinf, dMaxf;
    UInt64 segmentStart = 0;
    UInt64 segmentSize;

    if(nParticles == 0) {
        segmentSize = sizeof(UInt);
        memcpy(&nParticles, &buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;
    }

    segmentSize = sizeof(float);
    memcpy(&dMinf, &buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;
    memcpy(&dMaxf, &buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;

    segmentSize = nParticles * sizeof(UInt16);
    __NT_REQUIRE(segmentStart + segmentSize == buffer.size());
    StdVT_UInt16 compressedData(nParticles);
    memcpy(compressedData.data(), &buffer.data()[segmentStart], segmentSize);

    Real_t dMin = static_cast<Real_t>(dMinf);
    Real_t dMax = static_cast<Real_t>(dMaxf);
    decompress(dvec, dMin, dMax, compressedData);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void decompress(StdVT<StdVT<Real_t>>& dvec, const StdVT<Real_t>& dMin, const StdVT<Real_t>& dMax, const StdVT<StdVT_UInt16>& compressedData) {
    static_assert(std::is_floating_point_v<Real_t>);
    __NT_REQUIRE(compressedData.size() == dMin.size() && compressedData.size() == dMax.size());

    dvec.resize(compressedData.size());
    Scheduler::parallel_for(dvec.size(), [&](size_t i) { decompress(dvec[i], dMin[i], dMax[i], compressedData[i]); });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void decompress(StdVT<StdVT<Real_t>>& dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/) {
    static_assert(std::is_floating_point_v<Real_t>);
    StdVT<float> dMinf, dMaxf;

    UInt64 segmentStart = 0;
    UInt64 segmentSize;

    if(nParticles == 0) {
        segmentSize = sizeof(UInt);
        memcpy(&nParticles, &buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;
    }

    dMinf.resize(nParticles);
    dMaxf.resize(nParticles);
    StdVT<StdVT_UInt16> compressedData(nParticles);
    for(UInt i = 0; i < nParticles; ++i) {
        UInt iSize;
        segmentSize = sizeof(UInt);
        memcpy(&iSize, &buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;

        segmentSize = sizeof(float);
        memcpy(&dMinf[i], &buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;
        memcpy(&dMaxf[i], &buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;

        compressedData[i].resize(iSize);
        segmentSize = iSize * sizeof(UInt16);
        __NT_REQUIRE(segmentStart + segmentSize <= buffer.size());
        memcpy(compressedData[i].data(), &buffer.data()[segmentStart], segmentSize);
    }

    for(UInt i = 0; i < nParticles; ++i) {
        decompress(dvec[i], static_cast<Real_t>(dMinf[i]), static_cast<Real_t>(dMaxf[i]), compressedData[i]);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t> bool loadParticlesFromObj(const String& fileName, StdVT_VecX<N, Real_t>& positions) {
    std::ifstream file(fileName.c_str());
    if(!file.is_open()) {
        return false;
    }

    String line;
    while(std::getline(file, line)) {
        std::istringstream ls(line);
        String             label;
        ls >> label;
        if(label != "v") {
            continue;
        }
        if constexpr(N == 2) {
            Vec2f v;
            ls >> v.x >> v.y;
            positions.push_back(v);
        } else {
            Vec3f v;
            ls >> v.x >> v.y >> v.z;
            positions.push_back(v);
        }
    }
    file.close();
    return true;
}

template<Int N, class Real_t> bool loadParticlesFromBGEO(const String& fileName, StdVT_VecX<N, Real_t>& positions, Real_t& particleRadius) {
    Partio::ParticlesDataMutable* bgeoParticles = Partio::read(fileName.c_str());
    if(bgeoParticles == nullptr) { return false; }
    Partio::ParticleAttribute attrPosition;
    if(!bgeoParticles->attributeInfo("position", attrPosition)) {
        return false;
    }
    ////////////////////////////////////////////////////////////////////////////////
    if(Partio::ParticleAttribute attrRadius; bgeoParticles->attributeInfo("pscale", attrRadius)) {
        const float* radius = bgeoParticles->data<float>(attrRadius, 0);
        if(std::abs(particleRadius) < Real_t(1e-10)) {
            particleRadius = static_cast<Real_t>(radius[0]);
        } else {
            __NT_REQUIRE(std::abs(particleRadius - radius[0]) < MEpsilon<Real_t>());
        }
    }
    positions.resize(bgeoParticles->numParticles());
    auto* positionsPtr = reinterpret_cast<Real_t*>(positions.data());
    for(int p = 0; p < bgeoParticles->numParticles(); ++p) {
        const float* pos = bgeoParticles->data<float>(attrPosition, p);
        for(int i = 0; i < N; ++i) {
            positionsPtr[p * N + i] = static_cast<Real_t>(pos[i]);
        }
    }
    bgeoParticles->release();
    return true;
}

template<Int N, class Real_t> bool loadParticlesFromBNN(const String& fileName, StdVT_VecX<N, Real_t>& positions, Real_t& particleRadius) {
    return ParticleSerialization::loadParticle<N, Real_t>(fileName, positions, particleRadius);
}

template<Int N, class Real_t> bool loadParticlesFromBinary(const String& fileName, StdVT_VecX<N, Real_t>& positions, Real_t& particleRadius) {
    std::ifstream file(fileName.c_str(), std::ios::binary | std::ios::in);
    if(!file.is_open()) {
        return false;
    }

    UInt   nParticles;
    Real_t tmpRadius;
    file.read((char*)&nParticles, sizeof(UInt));
    file.read((char*)&tmpRadius,  sizeof(Real_t));
    if(particleRadius == 0) {
        particleRadius = tmpRadius;
    } else if(std::abs(particleRadius - tmpRadius) > MEpsilon<Real_t>()) {
        file.close();
        return false;
    }

    positions.resize(nParticles);
    file.read((char*)positions.data(), nParticles * sizeof(VecX<N, Real_t>));
    file.close();
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t> bool saveParticlesToObj(const String& fileName, const StdVT_VecX<N, Real_t>& positions) {
    std::ofstream file(fileName.c_str(), std::ios::out);
    if(!file.is_open()) {
        return false;
    }
    file << "# Num. particles: " << std::to_string(positions.size()) << "\n";
    const auto* positionsPtr = reinterpret_cast<const Real_t*>(positions.data());
    for(size_t p = 0; p < positions.size(); ++p) {
        file << "v";
        file << " " << std::to_string(positionsPtr[p * N]);
        file << " " << std::to_string(positionsPtr[p * N + 1]);
        if constexpr(N == 3) {
            file << " " << std::to_string(positionsPtr[p * N + 2]);
        } else {
            file << " 0";
        }
        file << "\n";
    }
    file.close();
    return true;
}

template<Int N, class Real_t> bool saveParticlesToBGEO(const String& fileName, const StdVT_VecX<N, Real_t>& positions, Real_t particleRadius) {
    Partio::ParticlesDataMutable* bgeoParticle = Partio::create();
    Partio::ParticleAttribute     attrRadius   = bgeoParticle->addAttribute("pscale", Partio::FLOAT, 1);
    Partio::ParticleAttribute     attrPosition = bgeoParticle->addAttribute("position", Partio::VECTOR, 3);

    const auto* positionsPtr = reinterpret_cast<const Real_t*>(positions.data());
    for(size_t p = 0; p < positions.size(); ++p) {
        size_t particle = bgeoParticle->addParticle();
        float* radius   = bgeoParticle->dataWrite<float>(attrRadius, particle);
        float* pos      = bgeoParticle->dataWrite<float>(attrPosition, particle);

        radius[0] = static_cast<float>(particleRadius);
        pos[0]    = static_cast<float>(positionsPtr[p * N]);
        pos[1]    = static_cast<float>(positionsPtr[p * N + 1]);
        if constexpr(N == 3) {
            pos[2] = static_cast<float>(positionsPtr[p * N + 2]);
        } else {
            pos[2] = 0;
        }
    }
    Partio::write(fileName.c_str(), *bgeoParticle);
    bgeoParticle->release();
    return true;
}

template<Int N, class Real_t> bool saveParticlesToBNN(const String& fileName, const StdVT_VecX<N, Real_t>& positions, Real_t particleRadius) {
    ParticleSerialization::saveParticle<N, Real_t>(fileName, positions, particleRadius);
    return true;
}

template<Int N, class Real_t> bool saveParticlesToBinary(const String& fileName, const StdVT_VecX<N, Real_t>& positions, Real_t particleRadius) {
    std::ofstream file(fileName.c_str(), std::ios::binary | std::ios::out);
    if(!file.is_open()) {
        return false;
    }

    UInt nParticles = static_cast<UInt>(positions.size());
    file.write((const char*)&nParticles,      sizeof(UInt));
    file.write((const char*)&particleRadius,  sizeof(Real_t));
    file.write((const char*)positions.data(), positions.size() * sizeof(VecX<N, Real_t>));
    file.close();
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// non-template functions
void connectedComponentAnalysis(const StdVT<StdVT_UInt>& connectionList, StdVT_Int8& componentIdx, UInt& nComponents) {
    componentIdx.assign(connectionList.size(), Int8(-1));

    // label from first particle
    UInt nProcessed = 0;
    nProcessed += spawnComponent(0, 0, 0, connectionList, componentIdx);

    nComponents = 1;
    Int8 currentCompIdx = 0;

    bool new_label = false;

    while(nProcessed < connectionList.size()) {
        bool bLabeled = false;

        for(UInt p = 0; p < connectionList.size(); ++p) {
            // Firstly, find any particle that has not been label
            // and that particle has a labeled neighbor
            if(componentIdx[p] >= 0) {
                continue;
            }

            if(new_label) {
                nProcessed += spawnComponent(p, 0, currentCompIdx, connectionList, componentIdx);

                bLabeled  = true;
                new_label = false;
            } else {
                // find the component index from neighbor
                Int8 neighborCompIdx = -1;

                for(UInt q : connectionList[p]) {
                    if(componentIdx[q] >= 0) {
                        neighborCompIdx = componentIdx[q];
                        break;
                    }
                }

                // has a labeled neighbor?
                // get component id from neighbor
                if(neighborCompIdx >= 0) {
                    nProcessed += spawnComponent(p, 0, neighborCompIdx, connectionList, componentIdx);
                    bLabeled    = true;
                }
            }
        }

        // not any particle has been labeled in the last loop
        // while num_process still < num particles
        // that means, we arrive at a new component
        if(!bLabeled) {
            ++currentCompIdx;
            ++nComponents;

            new_label = true;
        }
    }

    __NT_REQUIRE(nProcessed == connectionList.size());
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
UInt spawnComponent(UInt p, Int depth, UInt8 currentIdx, const StdVT<StdVT_UInt>& connectionList, StdVT_Int8& componentIdx) {
    componentIdx[p] = currentIdx;
    UInt nProcessed = 1;

    // max depth = 1024 to avoid stack overflow due to many time recursively call this function
    if(depth < 1024) {
        for(UInt q : connectionList[p]) {
            if(componentIdx[q] < 0) {
                nProcessed += spawnComponent(q, depth + 1, currentIdx, connectionList, componentIdx);
            }
        }
    }

    return nProcessed;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __BNN_INSTANTIATE_GET_AABB_COMMON_VEC_DIM(type)                                                     \
    template std::pair<VecX<2, type>, VecX<2, type>> getAABB<2, type>(const StdVT_VecX<2, type>&positions); \
    template std::pair<VecX<3, type>, VecX<3, type>> getAABB<3, type>(const StdVT_VecX<3, type>&positions); \
    template std::pair<VecX<4, type>, VecX<4, type>> getAABB<4, type>(const StdVT_VecX<4, type>&positions);

__BNN_INSTANTIATE_GET_AABB_COMMON_VEC_DIM(float)
__BNN_INSTANTIATE_GET_AABB_COMMON_VEC_DIM(double)

////////////////////////////////////////////////////////////////////////////////
#define __BNN_INSTANTIATE_GET_CENTER_COMMON_VEC_DIM(type)                           \
    template VecX<2, type> getCenter<2, type>(const StdVT_VecX<2, type>&positions); \
    template VecX<3, type> getCenter<3, type>(const StdVT_VecX<3, type>&positions); \
    template VecX<4, type> getCenter<4, type>(const StdVT_VecX<4, type>&positions);

__BNN_INSTANTIATE_GET_CENTER_COMMON_VEC_DIM(float)
__BNN_INSTANTIATE_GET_CENTER_COMMON_VEC_DIM(double)
////////////////////////////////////////////////////////////////////////////////
#define __BNN_INSTANTIATE_COMPRESS_COMMON_VEC_DIM(type)                                                                                     \
    template void compress<2, type>(const StdVT_VecX<2, type>&dvec, VecX<2, type>&dMin, VecX<2, type>&dMax, StdVT_UInt16 & compressedData); \
    template void compress<3, type>(const StdVT_VecX<3, type>&dvec, VecX<3, type>&dMin, VecX<3, type>&dMax, StdVT_UInt16 & compressedData); \
    template void compress<4, type>(const StdVT_VecX<4, type>&dvec, VecX<4, type>&dMin, VecX<4, type>&dMax, StdVT_UInt16 & compressedData); \
    template void compress<2, type>(const StdVT_VecX<2, type>&dvec, DataBuffer & buffer, bool bWriteVectorSize /*= true*/);                 \
    template void compress<3, type>(const StdVT_VecX<3, type>&dvec, DataBuffer & buffer, bool bWriteVectorSize /*= true*/);                 \
    template void compress<4, type>(const StdVT_VecX<4, type>&dvec, DataBuffer & buffer, bool bWriteVectorSize /*= true*/);                 \
    template void compress<2, type>(const StdVT<MatXxX<2, type>>&dvec, type & dMin, type & dMax, StdVT_UInt16 & compressedData);            \
    template void compress<3, type>(const StdVT<MatXxX<3, type>>&dvec, type & dMin, type & dMax, StdVT_UInt16 & compressedData);            \
    template void compress<4, type>(const StdVT<MatXxX<4, type>>&dvec, type & dMin, type & dMax, StdVT_UInt16 & compressedData);            \
    template void compress<2, type>(const StdVT<MatXxX<2, type>>&dvec, DataBuffer & buffer, bool bWriteVectorSize /*= true*/);              \
    template void compress<3, type>(const StdVT<MatXxX<3, type>>&dvec, DataBuffer & buffer, bool bWriteVectorSize /*= true*/);              \
    template void compress<4, type>(const StdVT<MatXxX<4, type>>&dvec, DataBuffer & buffer, bool bWriteVectorSize /*= true*/);

__BNN_INSTANTIATE_COMPRESS_COMMON_VEC_DIM(float)
__BNN_INSTANTIATE_COMPRESS_COMMON_VEC_DIM(double)
////////////////////////////////////////////////////////////////////////////////
#define __BNN_INSTANTIATE_COMPRESS(type)                                                                                                 \
    template void compress<type>(const StdVT<type>&dvec, DataBuffer & buffer, bool bWriteVectorSize /*= true*/);                         \
    template void compress<type>(const StdVT<type>&dvec, type & dMin, type & dMax, StdVT_UInt16 & compressedData);                       \
    template void compress<type>(const StdVT<StdVT<type>>&dvec, StdVT<type>&dMin, StdVT<type>&dMax, StdVT<StdVT_UInt16>&compressedData); \
    template void compress<type>(const StdVT<StdVT<type>>&dvec, DataBuffer & buffer, bool bWriteVectorSize /*= true*/);

__BNN_INSTANTIATE_COMPRESS(float)
__BNN_INSTANTIATE_COMPRESS(double)
////////////////////////////////////////////////////////////////////////////////
#define __BNN_INSTANTIATE_DECOMPRESS_COMMON_VEC_DIM(type)                                                                                                \
    template void decompress<2, type>(StdVT_VecX<2, type>&dvec, const VecX<2, type>&dMin, const VecX<2, type>&dMax, const StdVT_UInt16& compressedData); \
    template void decompress<3, type>(StdVT_VecX<3, type>&dvec, const VecX<3, type>&dMin, const VecX<3, type>&dMax, const StdVT_UInt16& compressedData); \
    template void decompress<4, type>(StdVT_VecX<4, type>&dvec, const VecX<4, type>&dMin, const VecX<4, type>&dMax, const StdVT_UInt16& compressedData); \
    template void decompress<2, type>(StdVT_VecX<2, type>&dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/);                                      \
    template void decompress<3, type>(StdVT_VecX<3, type>&dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/);                                      \
    template void decompress<4, type>(StdVT_VecX<4, type>&dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/);                                      \
    template void decompress<2, type>(StdVT<MatXxX<2, type>>&dvec, type dMin, type dMax, const StdVT_UInt16& compressedData);                            \
    template void decompress<3, type>(StdVT<MatXxX<3, type>>&dvec, type dMin, type dMax, const StdVT_UInt16& compressedData);                            \
    template void decompress<4, type>(StdVT<MatXxX<4, type>>&dvec, type dMin, type dMax, const StdVT_UInt16& compressedData);                            \
    template void decompress<2, type>(StdVT<MatXxX<2, type>>&dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/);                                   \
    template void decompress<3, type>(StdVT<MatXxX<3, type>>&dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/);                                   \
    template void decompress<4, type>(StdVT<MatXxX<4, type>>&dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/);

__BNN_INSTANTIATE_DECOMPRESS_COMMON_VEC_DIM(float)
__BNN_INSTANTIATE_DECOMPRESS_COMMON_VEC_DIM(double)
////////////////////////////////////////////////////////////////////////////////
#define __BNN_INSTANTIATE_DECOMPRESS(type)                                                                                                             \
    template void decompress<type>(StdVT<type>&dvec, type dMin, type dMax, const StdVT_UInt16& compressedData);                                        \
    template void decompress<type>(StdVT<type>&dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/);                                               \
    template void decompress<type>(StdVT<StdVT<type>>&dvec, const StdVT<type>&dMin, const StdVT<type>&dMax, const StdVT<StdVT_UInt16>&compressedData); \
    template void decompress<type>(StdVT<StdVT<type>>&dvec, const DataBuffer& buffer, UInt nParticles /*= 0*/);

__BNN_INSTANTIATE_DECOMPRESS(float)
__BNN_INSTANTIATE_DECOMPRESS(double)
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __BNN_INSTANTIATE_LOAD_SAVE_PARTICLES(N, type)                                                                            \
    template bool loadParticlesFromObj<N, type>(const String& fileName, StdVT_VecX<N, type>&positions);                           \
    template bool loadParticlesFromBGEO<N, type>(const String& fileName, StdVT_VecX<N, type>&positions, type & particleRadius);   \
    template bool loadParticlesFromBNN<N, type>(const String& fileName, StdVT_VecX<N, type>&positions, type & particleRadius);    \
    template bool loadParticlesFromBinary<N, type>(const String& fileName, StdVT_VecX<N, type>&positions, type & particleRadius); \
    template bool saveParticlesToObj<N, type>(const String& fileName, const StdVT_VecX<N, type>&positions);                       \
    template bool saveParticlesToBGEO<N, type>(const String& fileName, const StdVT_VecX<N, type>&positions, type particleRadius); \
    template bool saveParticlesToBNN<N, type>(const String& fileName, const StdVT_VecX<N, type>&positions, type particleRadius);  \
    template bool saveParticlesToBinary<N, type>(const String& fileName, const StdVT_VecX<N, type>&positions, type particleRadius);

__BNN_INSTANTIATE_LOAD_SAVE_PARTICLES(2, float)
__BNN_INSTANTIATE_LOAD_SAVE_PARTICLES(3, float)
__BNN_INSTANTIATE_LOAD_SAVE_PARTICLES(2, double)
__BNN_INSTANTIATE_LOAD_SAVE_PARTICLES(3, double)

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase::ParticleHelpers
