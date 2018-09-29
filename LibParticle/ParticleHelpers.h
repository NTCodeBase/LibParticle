//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//        __  __        _        __  ___ ____   __  ___
//       / / / /____ _ (_)_____ /  |/  // __ \ /  |/  /
//      / /_/ // __ `// // ___// /|_/ // /_/ // /|_/ /
//     / __  // /_/ // // /   / /  / // ____// /  / /
//    /_/ /_/ \__,_//_//_/   /_/  /_//_/    /_/  /_/
//
//    This file is part of HairMPM - Material Point Method for Hair Simulation.
//    Created: 2018. All rights reserved.
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#pragma once

#include <LibCommon/CommonSetup.h>
#include <LibCommon/Data/DataIO.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace ParticleHelpers
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// template functions are defined in a hpp file

template<Int N, class Real_t> std::pair<VecX<N, Real_t>, VecX<N, Real_t>> getAABB(const StdVT_VecX<N, Real_t>& positions);
template<Int N, class Real_t> VecX<N, Real_t>                               getCenter(const StdVT_VecX<N, Real_t>& positions);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t> void compress(const StdVT_VecX<N, Real_t>& dvec, VecX<N, Real_t>& dMin, VecX<N, Real_t>& dMax, StdVT_UInt16& compressedData);
template<Int N, class Real_t> void compress(const StdVT_VecX<N, Real_t>& dvec, DataBuffer& buffer, bool bWriteVectorSize = true);
template<Int N, class Real_t> void compress(const StdVT<MatXxX<N, Real_t>>& dvec, Real_t& dMin, Real_t& dMax, StdVT_UInt16& compressedData);
template<Int N, class Real_t> void compress(const StdVT<MatXxX<N, Real_t>>& dvec, DataBuffer& buffer, bool bWriteVectorSize = true);

template<class Real_t> void compress(const StdVT<Real_t>& dvec, Real_t& dMin, Real_t& dMax, StdVT_UInt16& compressedData);
template<class Real_t> void compress(const StdVT<Real_t>& dvec, DataBuffer& buffer, bool bWriteVectorSize = true);
template<class Real_t> void compress(const StdVT<StdVT<Real_t>>& dvec, StdVT<Real_t>& dMin, StdVT<Real_t>& dMax, StdVT<StdVT_UInt16>& compressedData);
template<class Real_t> void compress(const StdVT<StdVT<Real_t>>& dvec, DataBuffer& buffer, bool bWriteVectorSize = true);

template<Int N, class Real_t> void decompress(StdVT_VecX<N, Real_t>& dvec, const VecX<N, Real_t>& dMin, const VecX<N, Real_t>& dMax, const StdVT_UInt16& compressedData);
template<Int N, class Real_t> void decompress(StdVT_VecX<N, Real_t>& dvec, const DataBuffer& buffer, UInt nParticles = 0);
template<Int N, class Real_t> void decompress(StdVT<MatXxX<N, Real_t>>& dvec, Real_t dMin, Real_t dMax, const StdVT_UInt16& compressedData);
template<Int N, class Real_t> void decompress(StdVT<MatXxX<N, Real_t>>& dvec, const DataBuffer& buffer, UInt nParticles = 0);

template<class Real_t> void decompress(StdVT<Real_t>& dvec, Real_t dMin, Real_t dMax, const StdVT_UInt16& compressedData);
template<class Real_t> void decompress(StdVT<Real_t>& dvec, const DataBuffer& buffer, UInt nParticles = 0);
template<class Real_t> void decompress(StdVT<StdVT<Real_t>>& dvec, const StdVT<Real_t> dMin, const StdVT<Real_t>& dMax, const StdVT<StdVT_UInt16>& compressedData);
template<class Real_t> void decompress(StdVT<StdVT<Real_t>>& dvec, const DataBuffer& buffer, UInt nParticles = 0);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t> bool loadParticlesFromObj(const String& fileName, StdVT_VecX<N, Real_t>& positions);
template<Int N, class Real_t> bool loadParticlesFromBGEO(const String& fileName, StdVT_VecX<N, Real_t>& positions, Real_t& particleRadius);
template<Int N, class Real_t> bool loadParticlesFromBNN(const String& fileName, StdVT_VecX<N, Real_t>& positions, Real_t& particleRadius);
template<Int N, class Real_t> bool loadParticlesFromBinary(const String& fileName, StdVT_VecX<N, Real_t>& positions, Real_t& particleRadius);

template<Int N, class Real_t> bool saveParticlesToObj(const String& fileName, const StdVT_VecX<N, Real_t>& positions);
template<Int N, class Real_t> bool saveParticlesToBGEO(const String& fileName, const StdVT_VecX<N, Real_t>& positions, Real_t particleRadius);
template<Int N, class Real_t> bool saveParticlesToBNN(const String& fileName, const StdVT_VecX<N, Real_t>& positions, Real_t particleRadius);
template<Int N, class Real_t> bool saveParticlesToBinary(const String& fileName, const StdVT_VecX<N, Real_t>& positions, Real_t particleRadius);
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// non-template functions
void connectedComponentAnalysis(const StdVT<StdVT_UInt>& connectionList, StdVT_Int8& componentIdx, UInt& nComponents);
UInt spawnComponent(UInt p, Int depth, UInt8 currentIdx, const StdVT<StdVT_UInt>& connectionList, StdVT_Int8& componentIdx);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}   // end namespace ParticleHelpers
