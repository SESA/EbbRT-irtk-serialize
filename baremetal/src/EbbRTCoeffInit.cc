//          Copyright Boston University SESA Group 2013 - 2014.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
#include "EbbRTCoeffInit.h"
#include <vector>
#include <iostream>
#include <fstream>

#include <ebbrt/LocalIdMap.h>
#include <ebbrt/GlobalIdMap.h>
#include <ebbrt/UniqueIOBuf.h>
#include <ebbrt/IOBuf.h>
#include <ebbrt/Debug.h>
#include <ebbrt/Messenger.h>
#include <ebbrt/SpinBarrier.h>
#include <ebbrt/MulticoreEbb.h>
#include <ebbrt/EbbRef.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

EBBRT_PUBLISH_TYPE(, EbbRTCoeffInit);

using namespace ebbrt;

EbbRTCoeffInit& EbbRTCoeffInit::HandleFault(ebbrt::EbbId id) {
  {
    ebbrt::LocalIdMap::ConstAccessor accessor;
    auto found = ebbrt::local_id_map->Find(accessor, id);
    if (found) {
      auto& pr = *boost::any_cast<EbbRTCoeffInit*>(accessor->second);
      ebbrt::EbbRef<EbbRTCoeffInit>::CacheRef(id, pr);
      return pr;
    }
  }

  ebbrt::EventManager::EventContext context;
  auto f = ebbrt::global_id_map->Get(id);
  EbbRTCoeffInit* p;
  f.Then([&f, &context, &p, id](ebbrt::Future<std::string> inner) {
    p = new EbbRTCoeffInit(ebbrt::Messenger::NetworkId(inner.Get()), id);
    ebbrt::event_manager->ActivateContext(std::move(context));
  });
  ebbrt::event_manager->SaveContext(context);
  auto inserted = ebbrt::local_id_map->Insert(std::make_pair(id, p));
  if (inserted) {
    ebbrt::EbbRef<EbbRTCoeffInit>::CacheRef(id, *p);
    return *p;
  }

  delete p;
  // retry reading
  ebbrt::LocalIdMap::ConstAccessor accessor;
  ebbrt::local_id_map->Find(accessor, id);
  auto& pr = *boost::any_cast<EbbRTCoeffInit*>(accessor->second);
  ebbrt::EbbRef<EbbRTCoeffInit>::CacheRef(id, pr);
  return pr;
}

void EbbRTCoeffInit::Print(const char* str) {
  auto len = strlen(str) + 1;
  auto buf = ebbrt::MakeUniqueIOBuf(len);
  snprintf(reinterpret_cast<char*>(buf->MutData()), len, "%s", str);

  ebbrt::kprintf("Sending %d bytes\n", buf->ComputeChainDataLength());

  SendMessage(remote_nid_, std::move(buf));
}

void ResetOrigin(irtkGreyImage& image,
                 irtkRigidTransformation& transformation) {
  double ox, oy, oz;
  image.GetOrigin(ox, oy, oz);
  image.PutOrigin(0, 0, 0);
  transformation.PutTranslationX(ox);
  transformation.PutTranslationY(oy);
  transformation.PutTranslationZ(oz);
  transformation.PutRotationX(0);
  transformation.PutRotationY(0);
  transformation.PutRotationZ(0);
}

void ResetOrigin(irtkRealImage& image,
                 irtkRigidTransformation& transformation) {
  double ox, oy, oz;
  image.GetOrigin(ox, oy, oz);
  image.PutOrigin(0, 0, 0);
  transformation.PutTranslationX(ox);
  transformation.PutTranslationY(oy);
  transformation.PutTranslationZ(oz);
  transformation.PutRotationX(0);
  transformation.PutRotationY(0);
  transformation.PutRotationZ(0);
}

class MPMultiEbbCtr;
typedef EbbRef<MPMultiEbbCtr> MPMultiEbbCtrRef;

class MPMultiEbbCtrRoot {
 private:
  friend class MPMultiEbbCtr;
  EbbId _id;
  mutable MPMultiEbbCtr* repArray[ebbrt::Cpu::kMaxCpus];

  MPMultiEbbCtrRoot(EbbId id) : _id(id) { bzero(repArray, sizeof(repArray)); }

  void setRep(size_t i, MPMultiEbbCtr* rep) const { repArray[i] = rep; }

  int gatherVal(void) const;
  int otherGatherVal(void) const;
  void destroy() const;
};

class MPMultiEbbCtr : public MulticoreEbb<MPMultiEbbCtr, MPMultiEbbCtrRoot> {
  typedef MPMultiEbbCtrRoot Root;
  typedef MulticoreEbb<MPMultiEbbCtr, Root> Parent;
  const Root& _root;
  EbbId myId() { return _root._id; }

  int _val;
  MPMultiEbbCtr(const Root& root) : _root(root), _val(0) {
    _root.setRep(ebbrt::Cpu::GetMine(), this);
  }
  // give access to the constructor
  friend Parent;
  friend MPMultiEbbCtrRoot;

 public:
  void inc() { _val++; }
  void dec() { _val--; }
  void print() { ebbrt::kprintf("mycpu: %d\n", ebbrt::Cpu::GetMine()); }

  int val() { return _root.gatherVal(); }

  void destroy() { _root.destroy(); }

  static MPMultiEbbCtrRef Create(EbbId id = ebb_allocator->AllocateLocal()) {
    return Parent::Create(new Root(id), id);
  }
};

void MPMultiEbbCtrRoot::destroy(void) const {
  size_t numCores = ebbrt::Cpu::Count();
  for (size_t i = 0; numCores && i < ebbrt::Cpu::kMaxCpus; i++) {
    if (repArray[i]) {
      delete repArray[i];
      numCores--;
    }
  }
  delete this;
}

int MPMultiEbbCtrRoot::gatherVal(void) const {
  int gval = 0;
  size_t numCores = ebbrt::Cpu::Count();
  for (size_t i = 0; numCores && i < ebbrt::Cpu::kMaxCpus; i++) {
    if (repArray[i]) {
      gval = repArray[i]->_val;
      numCores--;
    }
  }
  return gval;
}

int MPMultiEbbCtrRoot::otherGatherVal(void) const {
  int gval = 0;
  LocalIdMap::ConstAccessor accessor;  // serves as a lock on the rep map
  auto found = local_id_map->Find(accessor, _id);
  if (!found)
    throw std::runtime_error("Failed to find root for MulticoreEbb");
  auto pair = boost::any_cast<std::pair<
      MPMultiEbbCtrRoot*, boost::container::flat_map<size_t, MPMultiEbbCtr*>>>(
      &accessor->second);
  const auto& rep_map = pair->second;
  for (auto it = rep_map.begin(); it != rep_map.end(); it++) {
    auto rep = boost::any_cast<const MPMultiEbbCtr*>(it->second);
    gval += rep->_val;
  }
  return gval;
};

void runCoeffInit(irtkRealImage& _mask, irtkRealImage& _reconstructed,
                  double& _quality_factor, size_t& _max_slices,
                  std::vector<irtkRealImage>& _slices,
                  std::vector<irtkRigidTransformation>& _transformations,
                  std::vector<SLICECOEFFS>& _volcoeffs,
                  std::vector<bool>& _slice_inside_cpu, int start, int end) {

  size_t inputIndex = 0;
  for (inputIndex = (size_t)start; inputIndex != (size_t)end; inputIndex++) {
    bool slice_inside;

    // get resolution of the volume
    double vx, vy, vz;
    _reconstructed.GetPixelSize(&vx, &vy, &vz);

    // volume is always isotropic
    double res = vx;

    // read the slice
    irtkRealImage& slice = _slices[inputIndex];

    // prepare structures for storage
    POINT3D p;
    VOXELCOEFFS empty;
    SLICECOEFFS slicecoeffs(slice.GetX(),
                            std::vector<VOXELCOEFFS>(slice.GetY(), empty));

    // to check whether the slice has an overlap with mask ROI
    slice_inside = false;

    // PSF will be calculated in slice space in higher resolution
    // get slice voxel size to define PSF
    double dx, dy, dz;
    slice.GetPixelSize(&dx, &dy, &dz);

    // sigma of 3D Gaussian (sinc with FWHM=dx or dy in-plane, Gaussian with
    // FWHM = dz through-plane)
    double sigmax = 1.2 * dx / 2.3548;
    double sigmay = 1.2 * dy / 2.3548;
    double sigmaz = dz / 2.3548;

    // calculate discretized PSF

    // isotropic voxel size of PSF - derived from resolution of reconstructed
    // volume
    double size = res / _quality_factor;

    // number of voxels in each direction
    // the ROI is 2*voxel dimension

    int xDim = round(2 * dx / size);
    int yDim = round(2 * dy / size);
    int zDim = round(2 * dz / size);

    // image corresponding to PSF
    irtkImageAttributes attr;
    attr._x = xDim;
    attr._y = yDim;
    attr._z = zDim;
    attr._dx = size;
    attr._dy = size;
    attr._dz = size;
    irtkRealImage PSF(attr);

    // centre of PSF
    double cx, cy, cz;
    cx = 0.5 * (xDim - 1);
    cy = 0.5 * (yDim - 1);
    cz = 0.5 * (zDim - 1);
    PSF.ImageToWorld(cx, cy, cz);

    double x, y, z;
    double sum = 0;
    int i, j, k;
    for (i = 0; i < xDim; i++)
      for (j = 0; j < yDim; j++)
        for (k = 0; k < zDim; k++) {
          x = i;
          y = j;
          z = k;
          PSF.ImageToWorld(x, y, z);
          x -= cx;
          y -= cy;
          z -= cz;
          // continuous PSF does not need to be normalized as discrete will be
          PSF(i, j, k) = exp(-x * x / (2 * sigmax * sigmax) -
                             y * y / (2 * sigmay * sigmay) -
                             z * z / (2 * sigmaz * sigmaz));
          sum += PSF(i, j, k);
        }
    PSF /= sum;

    // prepare storage for PSF transformed and resampled to the space of
    // reconstructed volume
    // maximum dim of rotated kernel - the next higher odd integer plus two to
    // accound for rounding error of tx,ty,tz.
    // Note conversion from PSF image coordinates to tPSF image coordinates
    // *size/res
    int dim =
        (floor(ceil(sqrt(double(xDim * xDim + yDim * yDim + zDim * zDim)) *
                    size / res) /
               2)) *
            2 +
        1 + 2;
    // prepare image attributes. Voxel dimension will be taken from the
    // reconstructed volume
    attr._x = dim;
    attr._y = dim;
    attr._z = dim;
    attr._dx = res;
    attr._dy = res;
    attr._dz = res;
    // create matrix from transformed PSF
    irtkRealImage tPSF(attr);
    // calculate centre of tPSF in image coordinates
    int centre = (dim - 1) / 2;

    // for each voxel in current slice calculate matrix coefficients
    int ii, jj, kk;
    int tx, ty, tz;
    int nx, ny, nz;
    int l, m, n;
    double weight;
    for (i = 0; i < slice.GetX(); i++)
      for (j = 0; j < slice.GetY(); j++)
        if (slice(i, j, 0) != -1) {
          // calculate centrepoint of slice voxel in volume space (tx,ty,tz)
          x = i;
          y = j;
          z = 0;
          slice.ImageToWorld(x, y, z);
          _transformations[inputIndex].Transform(x, y, z);
          _reconstructed.WorldToImage(x, y, z);
          tx = round(x);
          ty = round(y);
          tz = round(z);

          // Clear the transformed PSF
          for (ii = 0; ii < dim; ii++)
            for (jj = 0; jj < dim; jj++)
              for (kk = 0; kk < dim; kk++)
                tPSF(ii, jj, kk) = 0;

          // for each POINT3D of the PSF
          for (ii = 0; ii < xDim; ii++)
            for (jj = 0; jj < yDim; jj++)
              for (kk = 0; kk < zDim; kk++) {
                // Calculate the position of the POINT3D of
                // PSF centered over current slice voxel
                // This is a bit complicated because slices
                // can be oriented in any direction

                // PSF image coordinates
                x = ii;
                y = jj;
                z = kk;
                // change to PSF world coordinates - now real sizes in mm
                PSF.ImageToWorld(x, y, z);
                // centre around the centrepoint of the PSF
                x -= cx;
                y -= cy;
                z -= cz;

                // Need to convert (x,y,z) to slice image
                // coordinates because slices can have
                // transformations included in them (they are
                // nifti)  and those are not reflected in
                // PSF. In slice image coordinates we are
                // sure that z is through-plane

                // adjust according to voxel size
                x /= dx;
                y /= dy;
                z /= dz;
                // center over current voxel
                x += i;
                y += j;

                // convert from slice image coordinates to world coordinates
                slice.ImageToWorld(x, y, z);

                // x+=(vx-cx); y+=(vy-cy); z+=(vz-cz);
                // Transform to space of reconstructed volume
                _transformations[inputIndex].Transform(x, y, z);

                // Change to image coordinates
                _reconstructed.WorldToImage(x, y, z);

                // determine coefficients of volume voxels for position x,y,z
                // using linear interpolation

                // Find the 8 closest volume voxels

                // lowest corner of the cube
                nx = (int)floor(x);
                ny = (int)floor(y);
                nz = (int)floor(z);

                // not all neighbours might be in ROI, thus we need to
                // normalize
                //(l,m,n) are image coordinates of 8 neighbours in volume
                // space
                // for each we check whether it is in volume
                sum = 0;
                // to find wether the current slice voxel has overlap with ROI
                bool inside = false;
                for (l = nx; l <= nx + 1; l++)
                  if ((l >= 0) && (l < _reconstructed.GetX()))
                    for (m = ny; m <= ny + 1; m++)
                      if ((m >= 0) && (m < _reconstructed.GetY()))
                        for (n = nz; n <= nz + 1; n++)
                          if ((n >= 0) && (n < _reconstructed.GetZ())) {
                            weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) *
                                     (1 - fabs(n - z));
                            sum += weight;
                            if (_mask(l, m, n) == 1) {
                              inside = true;
                              slice_inside = true;
                            }
                          }
                // if there were no voxels do nothing
                if ((sum <= 0) || (!inside))
                  continue;
                // now calculate the transformed PSF
                for (l = nx; l <= nx + 1; l++)
                  if ((l >= 0) && (l < _reconstructed.GetX()))
                    for (m = ny; m <= ny + 1; m++)
                      if ((m >= 0) && (m < _reconstructed.GetY()))
                        for (n = nz; n <= nz + 1; n++)
                          if ((n >= 0) && (n < _reconstructed.GetZ())) {
                            weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) *
                                     (1 - fabs(n - z));

                            // image coordinates in tPSF
                            //(centre,centre,centre) in tPSF is aligned with
                            //(tx,ty,tz)
                            int aa, bb, cc;
                            aa = l - tx + centre;
                            bb = m - ty + centre;
                            cc = n - tz + centre;

                            // resulting value
                            double value = PSF(ii, jj, kk) * weight / sum;

                            // Check that we are in tPSF
                            if ((aa < 0) || (aa >= dim) || (bb < 0) ||
                                (bb >= dim) || (cc < 0) || (cc >= dim)) {
                              cerr << "Error while trying to populate tPSF. "
                                   << aa << " " << bb << " " << cc << endl;
                              cerr << l << " " << m << " " << n << endl;
                              cerr << tx << " " << ty << " " << tz << endl;
                              cerr << centre << endl;
                              tPSF.Write("tPSF.nii");
                              exit(1);
                            } else
                              // update transformed PSF
                              tPSF(aa, bb, cc) += value;
                          }
              }  // end of the loop for PSF points

          // store tPSF values
          for (ii = 0; ii < dim; ii++)
            for (jj = 0; jj < dim; jj++)
              for (kk = 0; kk < dim; kk++)
                if (tPSF(ii, jj, kk) > 0) {
                  p.x = ii + tx - centre;
                  p.y = jj + ty - centre;
                  p.z = kk + tz - centre;
                  p.value = tPSF(ii, jj, kk);
                  slicecoeffs[i][j].push_back(p);
                }
          // cout << " n = " << slicecoeffs[i][j].size() << std::endl;
        }  // end of loop for slice voxels

    _volcoeffs[inputIndex] = slicecoeffs;
    _slice_inside_cpu[inputIndex] = slice_inside;
  }
}

void EbbRTCoeffInit::doNothing() { ebbrt::kprintf("doNothing\n"); }

struct membuf : std::streambuf {
  membuf(char* begin, char* end) { this->setg(begin, begin, end); }
};

static size_t indexToCPU(size_t i) { return i; }

void EbbRTCoeffInit::ReceiveMessage(ebbrt::Messenger::NetworkId nid,
                                    std::unique_ptr<ebbrt::IOBuf>&& buffer) {
  auto output = std::string(reinterpret_cast<const char*>(buffer->Data()),
                            buffer->Length());
  // auto ip = nid.ToBytes().c_str();

  // ebbrt::kprintf("################# Serialization TEST
  // ###################\n");
  // ebbrt::kprintf("Received msg => %hhu:%hhu:%hhu:%hhu: %s\n",
  //               (unsigned char)ip[0], (unsigned char)ip[1],
  //               (unsigned char)ip[2], (unsigned char)ip[3], output.c_str());

  if (output[0] == 'G') {
    int start, end;
    // deserialize
    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};
    boost::archive::text_iarchive ia(stream);

    ia& start& end;

    std::vector<int> test(end - start);

    for (int k = 0; k < (end - start); k++) {
      ia& test[k];
    }

    ebbrt::kprintf("Received msg length: %d bytes\n", buffer->Length());
    ebbrt::kprintf("Number chain elements: %d\n", buffer->CountChainElements());
    ebbrt::kprintf("Computed chain length: %d bytes\n",
                   buffer->ComputeChainDataLength());

    // get number of cores/cpus on backend
    size_t ncpus = ebbrt::Cpu::Count();
    // get my cpu id
    size_t theCpu = ebbrt::Cpu::GetMine();

    // create a spin barrier on all cpus
    static ebbrt::SpinBarrier bar(ncpus);

    ebbrt::kprintf("number cpus: %d, mycpu: %d\n", ncpus, theCpu);

    // gets current context
    ebbrt::EventManager::EventContext context;

    // atomic type with value 0
    std::atomic<size_t> count(0);

    for (size_t i = 0; i < ncpus; i++) {
      // spawn jobs on each core using SpawnRemote
      ebbrt::event_manager->SpawnRemote(
          [theCpu, ncpus, &count, &test, i, &context, start, end]() {
            // get my cpu id
            size_t mycpu = ebbrt::Cpu::GetMine();

            int starte, ende, factor;
            factor = (int)ceil(test.size() / (float)ncpus);
            starte = i * factor;
            ende = i * factor + factor;
            ende = ((size_t)ende > test.size()) ? test.size() : ende;

            ebbrt::kprintf("theCpu: %d, mycpu: %d, start: %d, end: %d, "
                           "test.size(): %d, starte: %d, ende: %d\n",
                           theCpu, mycpu, start, end, test.size(), starte,
                           ende);

            for (int j = starte; j < ende; j++) {
              test[j] *= 2;
            }

            // atomically increment count
            count++;

            // barrier here ensures all cores run until this point
            bar.Wait();

            // basically wait until all cores reach this point
            while (count < (size_t)ncpus)
              ;

            // the cpu that initiated the SpawnRemote has the SaveContext
            if (mycpu == theCpu) {
              // activate context will return computation to instruction
              // after SaveContext below
              ebbrt::event_manager->ActivateContext(std::move(context));
            }
          },
          indexToCPU(
              i));  // if i don't add indexToCPU, one of the cores never run??
    }

    ebbrt::event_manager->SaveContext(context);

    ebbrt::kprintf("Context restored...\n");

    // serialize
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa& test& start& end;
    std::string ts = "G " + ofs.str();

    Print(ts.c_str());

  } else if (output[0] == 'Z') {
    irtkRealImage test;

    // deserialize
    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};
    boost::archive::text_iarchive ia(stream);
    ia& test;

    // serialize
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa& test;
    std::string ts = "E " + ofs.str();

    Print(ts.c_str());
  } else if (output[0] == 'D') {
    ebbrt::kprintf("Received msg length: %d bytes\n", buffer->Length());
    ebbrt::kprintf("Number chain elements: %d\n", buffer->CountChainElements());
    ebbrt::kprintf("Computed chain length: %d bytes\n",
                   buffer->ComputeChainDataLength());

    irtkRealImage test;

    // get number of cores/cpus on backend
    size_t ncpus = ebbrt::Cpu::Count();
    // get my cpu id
    size_t theCpu = ebbrt::Cpu::GetMine();

    // create a spin barrier on all cpus
    static ebbrt::SpinBarrier bar(ncpus);

    ebbrt::kprintf("number cpus: %d, mycpu: %d\n", ncpus, theCpu);

    // gets current context
    ebbrt::EventManager::EventContext context;

    // atomic type with value 0
    std::atomic<size_t> count(0);

    for (size_t i = 0; i < ncpus; i++) {
      // spawn jobs on each core using SpawnRemote
      ebbrt::event_manager->SpawnRemote(
          [theCpu, ncpus, &count, &context]() {
            // bar.Wait();
            // get my cpu id
            size_t mycpu = ebbrt::Cpu::GetMine();
            ebbrt::kprintf("theCpu: %d, mycpu: %d\n", theCpu, mycpu);
            // MPMultiEbbCtrRef ctr;
            // ctr = MPMultiEbbCtr::Create();
            // ebbrt::event_manager->Spawn([ctr]() {ctr->print();});
            // ctr->print();

            // atomically increment count
            count++;

            // barrier here ensures all cores run until this point
            bar.Wait();

            // basically wait until all cores reach this point
            while (count < (size_t)ncpus)
              ;

            // the cpu that initiated the SpawnRemote has the SaveContext
            if (ebbrt::Cpu::GetMine() == theCpu) {
              // activate context will return computation to instruction
              // after SaveContext below
              ebbrt::event_manager->ActivateContext(std::move(context));
            }
          },
          indexToCPU(
              i));  // if i don't add indexToCPU, one of the cores never run??
    }

    ebbrt::event_manager->SaveContext(context);

    ebbrt::kprintf("Context restored...\n");

    // deserialize
    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};
    boost::archive::text_iarchive ia(stream);
    ia& test;

    // serialize
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa& test;
    std::string ts = "D " + ofs.str();

    Print(ts.c_str());

  } else if (output[0] == 'A') {
    // deserialize and matrix mult
    std::istringstream ifs;
    ifs.str(output.substr(2, output.length() - 2));
    boost::archive::text_iarchive ia(ifs);
    irtkMatrix m1;
    irtkMatrix m2;
    irtkMatrix m3;
    ia& m1& m2;

    m3 = m1 * m2;

    // serialize m3
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa& m3;
    std::string ts = "A " + ofs.str();

    Print(ts.c_str());

    // ebbrt::kprintf("Received\n");

    // m1.Print();
    // m2.Print();

    // irtkMatrix m1;
    // std::vector<int> v2;
    // std::istringstream ifs;
    // ifs.str(output);
    // boost::archive::text_iarchive ia(ifs);
    // ia & v2;
    // for (auto &d: v2)
    // {
    //     ebbrt::kprintf("%d\n", d);
    // }

    // ebbrt::kprintf("Received vector, ready to shut down....\n");
  } else if (output[0] == 'B') {
    // deserialize and matrix mult
    std::istringstream ifs;
    ifs.str(output.substr(2, output.length() - 2));
    boost::archive::text_iarchive ia(ifs);
    irtkRealImage t;
    irtkRealImage _reconstructed;
    irtkImageAttributes attr;
    irtkRigidTransformation trans;
    ia& t& attr& _reconstructed& trans;

    irtkImageRigidRegistrationWithPadding registration;
    irtkGreyPixel smin, smax;
    irtkGreyImage target;
    irtkRealImage slice, w, b;
    double regCert;
    irtkResamplingWithPadding<irtkRealPixel> resampling(attr._dx, attr._dx,
                                                        attr._dx, -1);
    resampling.SetInput(&t);
    resampling.SetOutput(&t);
    resampling.Run();

    target = t;

    target.GetMinMax(&smin, &smax);

    if (smax > -1) {
      // put origin to zero
      irtkRigidTransformation offset;
      ResetOrigin(target, offset);

      irtkMatrix mo = offset.GetMatrix();
      irtkMatrix m = trans.GetMatrix();
      m = m * mo;
      trans.PutMatrix(m);

      irtkGreyImage source = _reconstructed;
      registration.SetInput(&target, &source);
      registration.SetOutput(&trans);
      registration.GuessParameterSliceToVolume();
      registration.SetTargetPadding(-1);
      registration.Run();

      regCert = registration.last_similarity;

      // undo the offset
      mo.Invert();
      m = trans.GetMatrix();
      m = m * mo;
      trans.PutMatrix(m);
    }

    // serialize m3
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa& t& trans& regCert;
    std::string ts = "B " + ofs.str();

    Print(ts.c_str());
  } else if (output[0] == 'C') {

    ebbrt::kprintf("Received msg length: %d bytes\n", buffer->Length());
    ebbrt::kprintf("Number chain elements: %d\n", buffer->CountChainElements());
    ebbrt::kprintf("Computed chain length: %d bytes\n",
                   buffer->ComputeChainDataLength());

    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};

    // deserialize
    // std::istringstream ifs;
    // std::string as(t);
    // ifs.str(as.substr(2, as.length() - 2));

    ebbrt::kprintf("Begin deserialization...\n");
    boost::archive::text_iarchive ia(stream);

    irtkRealImage _slice, _mask, _reconstructed;
    // irtkRigidTransformation _transformations;
    double _quality_factor;
    size_t _max_slices;
    size_t inputIndex = 0;
    std::vector<irtkRealImage> _slices;
    std::vector<irtkRigidTransformation> _transformations;
    std::vector<SLICECOEFFS> _volcoeffs;
    std::vector<bool> _slice_inside_cpu;

    // ia & _slice & _reconstructed & _quality_factor & _transformations &
    // _mask;

    ia& _reconstructed& _quality_factor& _mask& _max_slices& _slices&
    _transformations& _volcoeffs& _slice_inside_cpu;

    ebbrt::kprintf("Deserialized...\n");

    for (inputIndex = 0; inputIndex != _max_slices; inputIndex++) {
      bool slice_inside;

      // get resolution of the volume
      double vx, vy, vz;
      _reconstructed.GetPixelSize(&vx, &vy, &vz);

      // volume is always isotropic
      double res = vx;

      // read the slice
      irtkRealImage& slice = _slices[inputIndex];

      // prepare structures for storage
      POINT3D p;
      VOXELCOEFFS empty;
      SLICECOEFFS slicecoeffs(slice.GetX(),
                              std::vector<VOXELCOEFFS>(slice.GetY(), empty));

      // to check whether the slice has an overlap with mask ROI
      slice_inside = false;

      // PSF will be calculated in slice space in higher resolution
      // get slice voxel size to define PSF
      double dx, dy, dz;
      slice.GetPixelSize(&dx, &dy, &dz);

      // sigma of 3D Gaussian (sinc with FWHM=dx or dy in-plane, Gaussian with
      // FWHM = dz through-plane)
      double sigmax = 1.2 * dx / 2.3548;
      double sigmay = 1.2 * dy / 2.3548;
      double sigmaz = dz / 2.3548;

      // calculate discretized PSF

      // isotropic voxel size of PSF - derived from resolution of reconstructed
      // volume
      double size = res / _quality_factor;

      // number of voxels in each direction
      // the ROI is 2*voxel dimension

      int xDim = round(2 * dx / size);
      int yDim = round(2 * dy / size);
      int zDim = round(2 * dz / size);

      // image corresponding to PSF
      irtkImageAttributes attr;
      attr._x = xDim;
      attr._y = yDim;
      attr._z = zDim;
      attr._dx = size;
      attr._dy = size;
      attr._dz = size;
      irtkRealImage PSF(attr);

      // centre of PSF
      double cx, cy, cz;
      cx = 0.5 * (xDim - 1);
      cy = 0.5 * (yDim - 1);
      cz = 0.5 * (zDim - 1);
      PSF.ImageToWorld(cx, cy, cz);

      double x, y, z;
      double sum = 0;
      int i, j, k;
      for (i = 0; i < xDim; i++)
        for (j = 0; j < yDim; j++)
          for (k = 0; k < zDim; k++) {
            x = i;
            y = j;
            z = k;
            PSF.ImageToWorld(x, y, z);
            x -= cx;
            y -= cy;
            z -= cz;
            // continuous PSF does not need to be normalized as discrete will be
            PSF(i, j, k) = exp(-x * x / (2 * sigmax * sigmax) -
                               y * y / (2 * sigmay * sigmay) -
                               z * z / (2 * sigmaz * sigmaz));
            sum += PSF(i, j, k);
          }
      PSF /= sum;

      // prepare storage for PSF transformed and resampled to the space of
      // reconstructed volume
      // maximum dim of rotated kernel - the next higher odd integer plus two to
      // accound for rounding error of tx,ty,tz.
      // Note conversion from PSF image coordinates to tPSF image coordinates
      // *size/res
      int dim =
          (floor(ceil(sqrt(double(xDim * xDim + yDim * yDim + zDim * zDim)) *
                      size / res) /
                 2)) *
              2 +
          1 + 2;
      // prepare image attributes. Voxel dimension will be taken from the
      // reconstructed volume
      attr._x = dim;
      attr._y = dim;
      attr._z = dim;
      attr._dx = res;
      attr._dy = res;
      attr._dz = res;
      // create matrix from transformed PSF
      irtkRealImage tPSF(attr);
      // calculate centre of tPSF in image coordinates
      int centre = (dim - 1) / 2;

      // for each voxel in current slice calculate matrix coefficients
      int ii, jj, kk;
      int tx, ty, tz;
      int nx, ny, nz;
      int l, m, n;
      double weight;
      for (i = 0; i < slice.GetX(); i++)
        for (j = 0; j < slice.GetY(); j++)
          if (slice(i, j, 0) != -1) {
            // calculate centrepoint of slice voxel in volume space (tx,ty,tz)
            x = i;
            y = j;
            z = 0;
            slice.ImageToWorld(x, y, z);
            _transformations[inputIndex].Transform(x, y, z);
            _reconstructed.WorldToImage(x, y, z);
            tx = round(x);
            ty = round(y);
            tz = round(z);

            // Clear the transformed PSF
            for (ii = 0; ii < dim; ii++)
              for (jj = 0; jj < dim; jj++)
                for (kk = 0; kk < dim; kk++)
                  tPSF(ii, jj, kk) = 0;

            // for each POINT3D of the PSF
            for (ii = 0; ii < xDim; ii++)
              for (jj = 0; jj < yDim; jj++)
                for (kk = 0; kk < zDim; kk++) {
                  // Calculate the position of the POINT3D of
                  // PSF centered over current slice voxel
                  // This is a bit complicated because slices
                  // can be oriented in any direction

                  // PSF image coordinates
                  x = ii;
                  y = jj;
                  z = kk;
                  // change to PSF world coordinates - now real sizes in mm
                  PSF.ImageToWorld(x, y, z);
                  // centre around the centrepoint of the PSF
                  x -= cx;
                  y -= cy;
                  z -= cz;

                  // Need to convert (x,y,z) to slice image
                  // coordinates because slices can have
                  // transformations included in them (they are
                  // nifti)  and those are not reflected in
                  // PSF. In slice image coordinates we are
                  // sure that z is through-plane

                  // adjust according to voxel size
                  x /= dx;
                  y /= dy;
                  z /= dz;
                  // center over current voxel
                  x += i;
                  y += j;

                  // convert from slice image coordinates to world coordinates
                  slice.ImageToWorld(x, y, z);

                  // x+=(vx-cx); y+=(vy-cy); z+=(vz-cz);
                  // Transform to space of reconstructed volume
                  _transformations[inputIndex].Transform(x, y, z);

                  // Change to image coordinates
                  _reconstructed.WorldToImage(x, y, z);

                  // determine coefficients of volume voxels for position x,y,z
                  // using linear interpolation

                  // Find the 8 closest volume voxels

                  // lowest corner of the cube
                  nx = (int)floor(x);
                  ny = (int)floor(y);
                  nz = (int)floor(z);

                  // not all neighbours might be in ROI, thus we need to
                  // normalize
                  //(l,m,n) are image coordinates of 8 neighbours in volume
                  // space
                  // for each we check whether it is in volume
                  sum = 0;
                  // to find wether the current slice voxel has overlap with ROI
                  bool inside = false;
                  for (l = nx; l <= nx + 1; l++)
                    if ((l >= 0) && (l < _reconstructed.GetX()))
                      for (m = ny; m <= ny + 1; m++)
                        if ((m >= 0) && (m < _reconstructed.GetY()))
                          for (n = nz; n <= nz + 1; n++)
                            if ((n >= 0) && (n < _reconstructed.GetZ())) {
                              weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) *
                                       (1 - fabs(n - z));
                              sum += weight;
                              if (_mask(l, m, n) == 1) {
                                inside = true;
                                slice_inside = true;
                              }
                            }
                  // if there were no voxels do nothing
                  if ((sum <= 0) || (!inside))
                    continue;
                  // now calculate the transformed PSF
                  for (l = nx; l <= nx + 1; l++)
                    if ((l >= 0) && (l < _reconstructed.GetX()))
                      for (m = ny; m <= ny + 1; m++)
                        if ((m >= 0) && (m < _reconstructed.GetY()))
                          for (n = nz; n <= nz + 1; n++)
                            if ((n >= 0) && (n < _reconstructed.GetZ())) {
                              weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) *
                                       (1 - fabs(n - z));

                              // image coordinates in tPSF
                              //(centre,centre,centre) in tPSF is aligned with
                              //(tx,ty,tz)
                              int aa, bb, cc;
                              aa = l - tx + centre;
                              bb = m - ty + centre;
                              cc = n - tz + centre;

                              // resulting value
                              double value = PSF(ii, jj, kk) * weight / sum;

                              // Check that we are in tPSF
                              if ((aa < 0) || (aa >= dim) || (bb < 0) ||
                                  (bb >= dim) || (cc < 0) || (cc >= dim)) {
                                cerr << "Error while trying to populate tPSF. "
                                     << aa << " " << bb << " " << cc << endl;
                                cerr << l << " " << m << " " << n << endl;
                                cerr << tx << " " << ty << " " << tz << endl;
                                cerr << centre << endl;
                                tPSF.Write("tPSF.nii");
                                exit(1);
                              } else
                                // update transformed PSF
                                tPSF(aa, bb, cc) += value;
                            }
                }  // end of the loop for PSF points

            // store tPSF values
            for (ii = 0; ii < dim; ii++)
              for (jj = 0; jj < dim; jj++)
                for (kk = 0; kk < dim; kk++)
                  if (tPSF(ii, jj, kk) > 0) {
                    p.x = ii + tx - centre;
                    p.y = jj + ty - centre;
                    p.z = kk + tz - centre;
                    p.value = tPSF(ii, jj, kk);
                    slicecoeffs[i][j].push_back(p);
                  }
            // cout << " n = " << slicecoeffs[i][j].size() << std::endl;
          }  // end of loop for slice voxels

      _volcoeffs[inputIndex] = slicecoeffs;
      _slice_inside_cpu[inputIndex] = slice_inside;
    }

    // serialize m3
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);

    oa& _slices& _transformations& _volcoeffs& _slice_inside_cpu;
    // oa& _slice& _transformations& slicecoeffs& slice_inside;

    std::string ts = "C " + ofs.str();
    ebbrt::kprintf("Sending length: %d\n", ts.length());
    Print(ts.c_str());
  } else if (output[0] == 'F') {
    ebbrt::kprintf("Received msg length: %d bytes\n", buffer->Length());
    ebbrt::kprintf("Number chain elements: %d\n", buffer->CountChainElements());
    ebbrt::kprintf("Computed chain length: %d bytes\n",
                   buffer->ComputeChainDataLength());

    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};

    ebbrt::kprintf("Begin deserialization...\n");
    boost::archive::text_iarchive ia(stream);

    irtkRealImage _slice, _mask, _reconstructed;
    double _quality_factor;
    size_t _max_slices;
    // size_t inputIndex = 0;
    std::vector<irtkRealImage> _slices;
    std::vector<irtkRigidTransformation> _transformations;
    std::vector<SLICECOEFFS> _volcoeffs;
    std::vector<bool> _slice_inside_cpu;

    ia& _reconstructed& _quality_factor& _mask& _max_slices& _slices&
    _transformations& _volcoeffs& _slice_inside_cpu;

    ebbrt::kprintf("Deserialized...\n");

    // get number of cores/cpus on backend
    size_t ncpus = ebbrt::Cpu::Count();

    // create a spin barrier on all cpus
    static ebbrt::SpinBarrier bar(ncpus);

    ebbrt::kprintf("number cpu: %d\n", ncpus);

    // gets current context
    ebbrt::EventManager::EventContext context;

    // atomic type with value 0
    std::atomic<size_t> count(0);

    // get my cpu id
    //    size_t theCpu = ebbrt::Cpu::GetMine();

    /*
        for (size_t i = 0; i < ncpus; i++) {
          // spawn jobs on each core using SpawnRemote
          ebbrt::event_manager->SpawnRemote(
              [theCpu, ncpus, &count, &context, &_reconstructed,
       &_quality_factor,
               &_mask, &_max_slices, &_slices, &_transformations, &_volcoeffs,
               &_slice_inside_cpu, i]() {
                // get my cpu id
                size_t mycpu = ebbrt::Cpu::GetMine();

                int start, end, factor;
                factor = _max_slices / (int)ncpus;
                start = i * factor;
                end = i * factor + factor;
                end = ((size_t)end > _max_slices) ? _max_slices : end;

                ebbrt::kprintf("theCpu: %d, mycpu: %d, start: %d, end: %d\n",
                               theCpu, mycpu, start, end);

                //runCoeffInit(_mask, _reconstructed, _quality_factor,
       _max_slices,
                //           _slices, _transformations, _volcoeffs,
                //           _slice_inside_cpu, start, end);

                // atomically increment count
                count++;

                // barrier here ensures all cores run until this point
                bar.Wait();

                // basically wait until all cores reach this point
                while (count < (size_t)ncpus)
                  ;

                // the cpu that initiated the SpawnRemote has the SaveContext
                if (mycpu == theCpu) {
                  // activate context will return computation to instruction
                  // after SaveContext below
                  ebbrt::event_manager->ActivateContext(std::move(context));
                }
              },
              indexToCPU(
                  i));  // if i don't add indexToCPU, one of the cores never
       run??
        }*/

    ebbrt::event_manager->SaveContext(context);

    ebbrt::kprintf("Context restored...\n");
    // serialize m3
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);

    oa& _slices& _transformations& _volcoeffs& _slice_inside_cpu;
    // oa& _slice& _transformations& slicecoeffs& slice_inside;

    std::string ts = "C " + ofs.str();
    Print(ts.c_str());
  } else if (output[0] == 'E') {
    ebbrt::kprintf("Received msg length: %d bytes\n", buffer->Length());
    ebbrt::kprintf("Number chain elements: %d\n", buffer->CountChainElements());
    ebbrt::kprintf("Computed chain length: %d bytes\n",
                   buffer->ComputeChainDataLength());

    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};

    ebbrt::kprintf("Begin deserialization...\n");
    boost::archive::text_iarchive ia(stream);

    int start, end, diff;
    ia& start& end;
    diff = end - start;

    std::vector<irtkRealImage> _slices(diff);
    std::vector<irtkRigidTransformation> _transformations(diff);
    std::vector<SLICECOEFFS> _volcoeffs(diff);
    std::vector<bool> _slice_inside_cpu;

    irtkRealImage _slice, _mask, _reconstructed;
    double _quality_factor;
    size_t _max_slices;

    for (int k = 0; k < diff; k++) {
      ia& _slices[k];
    }

    for (int k = 0; k < diff; k++) {
      ia& _transformations[k];
    }

    for (int k = 0; k < diff; k++) {
      ia& _volcoeffs[k];
    }

    // ia & _slices & _transformations & _volcoeffs & _slice_inside_cpu;

    ia& _slice_inside_cpu;

    //    for (int k = 0; k < diff; k++) {
    //    ia& _slice_inside_cpu[k];
    //}

    ia& _reconstructed& _quality_factor& _mask& _max_slices;

    ebbrt::kprintf("Deserialized...\n");

    ebbrt::kprintf("_slices: %d\n", _slices.size());
    ebbrt::kprintf("_transformations: %d\n", _transformations.size());
    ebbrt::kprintf("_volcoeffs: %d\n", _volcoeffs.size());
    ebbrt::kprintf("_slice_inside_cpu: %d\n", _slice_inside_cpu.size());

    // int ende = (int)_slices.size();

    // runCoeffInit(_mask, _reconstructed, _quality_factor, _max_slices,
    //		 _slices, _transformations, _volcoeffs,
    //		 _slice_inside_cpu, 0, diff);

    // get number of cores/cpus on backend
    size_t ncpus = ebbrt::Cpu::Count();

    // create a spin barrier on all cpus
    static ebbrt::SpinBarrier bar(ncpus);

    ebbrt::kprintf("number cpu: %d\n", ncpus);

    // gets current context
    ebbrt::EventManager::EventContext context;

    // atomic type with value 0
    std::atomic<size_t> count(0);

    // get my cpu id
    size_t theCpu = ebbrt::Cpu::GetMine();

    for (size_t i = 0; i < ncpus; i++) {
      // spawn jobs on each core using SpawnRemote
      ebbrt::event_manager->SpawnRemote(
          [theCpu, ncpus, &count, &context, &_reconstructed, &_quality_factor,
           &_mask, &_max_slices, &_slices, &_transformations, &_volcoeffs,
           &_slice_inside_cpu, i, diff]() {
            // get my cpu id
            size_t mycpu = ebbrt::Cpu::GetMine();

            int starte, ende, factor;
            factor = (int)ceil(diff / (float)ncpus);
            starte = i * factor;
            ende = i * factor + factor;
            ende = (ende > diff) ? diff : ende;

            // ebbrt::kprintf("theCpu: %d, mycpu: %d, start: %d, end:, theCpu,
            // mycpu, start, end);

            runCoeffInit(_mask, _reconstructed, _quality_factor, _max_slices,
                         _slices, _transformations, _volcoeffs,
                         _slice_inside_cpu, starte, ende);

            // atomically increment count
            count++;

            // barrier here ensures all cores run until this point
            bar.Wait();

            // basically wait until all cores reach this point
            while (count < (size_t)ncpus)
              ;

            // the cpu that initiated the SpawnRemote has the SaveContext
            if (mycpu == theCpu) {
              // activate context will return computation to instruction
              // after SaveContext below
              ebbrt::event_manager->ActivateContext(std::move(context));
            }
          },
          indexToCPU(
              i));  // if i don't add indexToCPU, one of the cores never run ? ?
    }

    ebbrt::event_manager->SaveContext(context);

    ebbrt::kprintf("Context restored...\n");

    // serialize m3
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);

    ebbrt::kprintf("start: %d end: %d\n", start, end);
    ebbrt::kprintf("_slices: %d\n", _slices.size());
    ebbrt::kprintf("_transformations: %d\n", _transformations.size());
    ebbrt::kprintf("_volcoeffs: %d\n", _volcoeffs.size());
    ebbrt::kprintf("_slice_inside_cpu: %d\n", _slice_inside_cpu.size());

    oa& start& end;

    for (int j = 0; j < diff; j++) {
      oa& _slices[j];
    }

    for (int j = 0; j < diff; j++) {
      oa& _transformations[j];
    }

    for (int j = 0; j < diff; j++) {
      oa& _volcoeffs[j];
    }

    //& _slices& _transformations& _volcoeffs& _slice_inside_cpu;
    oa& _slice_inside_cpu;

    std::string ts = "E " + ofs.str();
    ebbrt::kprintf("ts length: %d\n", ts.length());

    Print(ts.c_str());
  }
}
