//          Copyright Boston University SESA Group 2013 - 2014.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
#include "Printer.h"
#include <vector>
#include <iostream>
#include <fstream>

#include <ebbrt/LocalIdMap.h>
#include <ebbrt/GlobalIdMap.h>
#include <ebbrt/UniqueIOBuf.h>
#include <ebbrt/Debug.h>
#include <ebbrt/Messenger.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

EBBRT_PUBLISH_TYPE(, Printer);

Printer::Printer(ebbrt::Messenger::NetworkId nid)
    : Messagable<Printer>(kPrinterEbbId), remote_nid_(std::move(nid)) {}

Printer& Printer::HandleFault(ebbrt::EbbId id) {
  {
    ebbrt::LocalIdMap::ConstAccessor accessor;
    auto found = ebbrt::local_id_map->Find(accessor, id);
    if (found) {
      auto& pr = *boost::any_cast<Printer*>(accessor->second);
      ebbrt::EbbRef<Printer>::CacheRef(id, pr);
      return pr;
    }
  }

  ebbrt::EventManager::EventContext context;
  auto f = ebbrt::global_id_map->Get(id);
  Printer* p;
  f.Then([&f, &context, &p](ebbrt::Future<std::string> inner) {
    p = new Printer(ebbrt::Messenger::NetworkId(inner.Get()));
    ebbrt::event_manager->ActivateContext(std::move(context));
  });
  ebbrt::event_manager->SaveContext(context);
  auto inserted = ebbrt::local_id_map->Insert(std::make_pair(id, p));
  if (inserted) {
    ebbrt::EbbRef<Printer>::CacheRef(id, *p);
    return *p;
  }

  delete p;
  // retry reading
  ebbrt::LocalIdMap::ConstAccessor accessor;
  ebbrt::local_id_map->Find(accessor, id);
  auto& pr = *boost::any_cast<Printer*>(accessor->second);
  ebbrt::EbbRef<Printer>::CacheRef(id, pr);
  return pr;
}

void Printer::Print(const char* str) {
  auto len = strlen(str) + 1;
  auto buf = ebbrt::MakeUniqueIOBuf(len);
  snprintf(reinterpret_cast<char*>(buf->MutData()), len, "%s", str);
  SendMessage(remote_nid_, std::move(buf));
}

void ResetOrigin(irtkGreyImage &image,
		 irtkRigidTransformation &transformation) {
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

void ResetOrigin(irtkRealImage &image,
		 irtkRigidTransformation &transformation) {
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

void Printer::ReceiveMessage(ebbrt::Messenger::NetworkId nid,
                             std::unique_ptr<ebbrt::IOBuf>&& buffer) {
  auto output = std::string(reinterpret_cast<const char*>(buffer->Data()),
                            buffer->Length());
  auto ip = nid.ToBytes().c_str();

  ebbrt::kprintf("################# Serialization TEST ###################\n");
  ebbrt::kprintf("Received msg => %hhu:%hhu:%hhu:%hhu: %s\n",
                 (unsigned char)ip[0], (unsigned char)ip[1],
                 (unsigned char)ip[2], (unsigned char)ip[3], output.c_str());

  if (output[0] == 'A') {
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
    oa & t & trans & regCert;
    std::string ts = "B " + ofs.str();

    Print(ts.c_str());
  }
  else if (output[0] == 'C') {
      
      // deserialize
      std::istringstream ifs;
      ifs.str(output.substr(2, output.length() - 2));
      boost::archive::text_iarchive ia(ifs);

      irtkRealImage _slice;//, _mask, _reconstructed;
//      irtkRigidTransformation _transformations;
      //    double _quality_factor;
      ia & _slice;// & _reconstructed & _quality_factor & _transformations & _mask;
      
      //bool slice_inside;
      
      // get resolution of the volume
      //double vx, vy, vz;
      //_reconstructed.GetPixelSize(&vx, &vy, &vz);
      
      // volume is always isotropic
      //double res = vx;

      // read the slice
      //irtkRealImage &slice = _slice;

      // prepare structures for storage
      //POINT3D p;
      //VOXELCOEFFS empty;
      //SLICECOEFFS slicecoeffs(slice.GetX(),
      //                      std::vector<VOXELCOEFFS>(slice.GetY(), empty));

      // to check whether the slice has an overlap with mask ROI
      //slice_inside = false;
      /*
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
            _transformations.Transform(x, y, z);
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
                  _transformations.Transform(x, y, z);

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
                        if ((m >= 0) &&
                            (m < _reconstructed.GetY()))
                          for (n = nz; n <= nz + 1; n++)
                            if ((n >= 0) &&
                                (n < _reconstructed.GetZ())) {
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
                        if ((m >= 0) &&
                            (m < _reconstructed.GetY()))
                          for (n = nz; n <= nz + 1; n++)
                            if ((n >= 0) &&
                                (n < _reconstructed.GetZ())) {
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
                } // end of the loop for PSF points

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
          } // end of loop for slice voxels

      //reconstructor->_volcoeffs[inputIndex] = slicecoeffs;
      //reconstructor->_slice_inside_cpu[inputIndex] = slice_inside;
      */
      // serialize m3
      std::ostringstream ofs;
      boost::archive::text_oarchive oa(ofs);
      //oa & _slice & _transformations & slicecoeffs & slice_inside;
      oa & _slice;//_reconstructed;
      std::string ts = "C " + ofs.str();
      Print(ts.c_str());
  }
}
