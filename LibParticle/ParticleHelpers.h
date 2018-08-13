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

#include <CommonSetup.h>
#include <Data/DataIO.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace ParticleHelpers
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// template functions are defined in a hpp file

template<Int N, class RealType> std::pair<VecX<N, RealType>, VecX<N, RealType>> getAABB(const StdVT_VecX<N, RealType>& positions);
template<Int N, class RealType> VecX<N, RealType>                               getCenter(const StdVT_VecX<N, RealType>& positions);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType> void compress(const StdVT_VecX<N, RealType>& dvec, VecX<N, RealType>& dMin, VecX<N, RealType>& dMax, StdVT_UInt16& compressedData);
template<Int N, class RealType> void compress(const StdVT_VecX<N, RealType>& dvec, DataBuffer& buffer, bool bWriteVectorSize = true);
template<Int N, class RealType> void compress(const StdVT<MatXxX<N, RealType>>& dvec, RealType& dMin, RealType& dMax, StdVT_UInt16& compressedData);
template<Int N, class RealType> void compress(const StdVT<MatXxX<N, RealType>>& dvec, DataBuffer& buffer, bool bWriteVectorSize = true);

template<class RealType> void compress(const StdVT<RealType>& dvec, RealType& dMin, RealType& dMax, StdVT_UInt16& compressedData);
template<class RealType> void compress(const StdVT<RealType>& dvec, DataBuffer& buffer, bool bWriteVectorSize = true);
template<class RealType> void compress(const StdVT<StdVT<RealType>>& dvec, StdVT<RealType>& dMin, StdVT<RealType>& dMax, StdVT<StdVT_UInt16>& compressedData);
template<class RealType> void compress(const StdVT<StdVT<RealType>>& dvec, DataBuffer& buffer, bool bWriteVectorSize = true);

template<Int N, class RealType> void decompress(StdVT_VecX<N, RealType>& dvec, const VecX<N, RealType>& dMin, const VecX<N, RealType>& dMax, const StdVT_UInt16& compressedData);
template<Int N, class RealType> void decompress(StdVT_VecX<N, RealType>& dvec, const DataBuffer& buffer, UInt nParticles = 0);
template<Int N, class RealType> void decompress(StdVT<MatXxX<N, RealType>>& dvec, RealType dMin, RealType dMax, const StdVT_UInt16& compressedData);
template<Int N, class RealType> void decompress(StdVT<MatXxX<N, RealType>>& dvec, const DataBuffer& buffer, UInt nParticles = 0);

template<class RealType> void decompress(StdVT<RealType>& dvec, RealType dMin, RealType dMax, const StdVT_UInt16& compressedData);
template<class RealType> void decompress(StdVT<RealType>& dvec, const DataBuffer& buffer, UInt nParticles = 0);
template<class RealType> void decompress(StdVT<StdVT<RealType>>& dvec, const StdVT<RealType> dMin, const StdVT<RealType>& dMax, const StdVT<StdVT_UInt16>& compressedData);
template<class RealType> void decompress(StdVT<StdVT<RealType>>& dvec, const DataBuffer& buffer, UInt nParticles = 0);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType> bool loadParticlesFromObj(const String& fileName, StdVT_VecX<N, RealType>& positions);
template<Int N, class RealType> bool loadParticlesFromBGEO(const String& fileName, StdVT_VecX<N, RealType>& positions, RealType& particleRadius);
template<Int N, class RealType> bool loadParticlesFromBNN(const String& fileName, StdVT_VecX<N, RealType>& positions, RealType& particleRadius);
template<Int N, class RealType> bool loadParticlesFromBinary(const String& fileName, StdVT_VecX<N, RealType>& positions, RealType& particleRadius);

template<Int N, class RealType> bool saveParticlesToObj(const String& fileName, const StdVT_VecX<N, RealType>& positions);
template<Int N, class RealType> bool saveParticlesToBGEO(const String& fileName, const StdVT_VecX<N, RealType>& positions, RealType particleRadius);
template<Int N, class RealType> bool saveParticlesToBNN(const String& fileName, const StdVT_VecX<N, RealType>& positions, RealType particleRadius);
template<Int N, class RealType> bool saveParticlesToBinary(const String& fileName, const StdVT_VecX<N, RealType>& positions, RealType particleRadius);
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// non-template functions
void connectedComponentAnalysis(const StdVT<StdVT_UInt>& connectionList, StdVT_Int8& componentIdx, UInt& nComponents);
UInt spawnComponent(UInt p, Int depth, UInt8 currentIdx, const StdVT<StdVT_UInt>& connectionList, StdVT_Int8& componentIdx);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}   // end namespace ParticleHelpers
