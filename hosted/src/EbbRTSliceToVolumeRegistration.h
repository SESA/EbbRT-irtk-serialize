#ifndef _EBBRT_SLICE_TO_VOLUME_REGISTRATION_H_
#define _EBBRT_SLICE_TO_VOLUME_REGISTRATION_H_

#define NOMINMAX
#define _USE_MATH_DEFINES

#include <irtkReconstructionGPU.h>
#include <irtkResampling.h>
#include <irtkRegistration.h>
#include <irtkImageRigidRegistration.h>
#include <irtkImageRigidRegistrationWithPadding.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>
#include <math.h>
#include <stdlib.h>
#include <utils.h>
#include <string>

#include <ebbrt/GlobalIdMap.h>
#include <ebbrt/StaticIds.h>
#include <ebbrt/EbbRef.h>
#include <ebbrt/SharedEbb.h>
#include <ebbrt/Message.h>
#include <ebbrt/UniqueIOBuf.h>
#include "../../src/StaticEbbIds.h"

#include <iostream>

#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using namespace ebbrt;

class EbbRTSliceToVolumeRegistration;
typedef EbbRef<EbbRTSliceToVolumeRegistration>
EbbRTSliceToVolumeRegistrationEbbRef;

class EbbRTSliceToVolumeRegistration
    : public SharedEbb<EbbRTSliceToVolumeRegistration>,
      public Messagable<EbbRTSliceToVolumeRegistration> {

  irtkReconstruction* reconstructor;
  EbbId ebbid;
  std::string bmnid;
  irtkMatrix mym;
  irtkRealImage t;
  double regCert;
  irtkRigidTransformation trans;
  
  // this is used to save and load context
  ebbrt::EventManager::EventContext* emec{nullptr};

  // a void promise where we want to start some computation after some other
  // function has finished
  ebbrt::Promise<void> mypromise;

 public:
  EbbRTSliceToVolumeRegistration(irtkReconstruction* _reconstructor,
                                 EbbId _ebbid)
      : Messagable<EbbRTSliceToVolumeRegistration>(_ebbid),
        reconstructor(_reconstructor) {
    ebbid = _ebbid;
    bmnid = "";
  }

  // Create function that returns a future ebbref
  static ebbrt::Future<EbbRTSliceToVolumeRegistrationEbbRef>
  Create(irtkReconstruction* _reconstructor) {
    // calls the constructor in this class and returns a EbbRef
    auto ebbref = SharedEbb<EbbRTSliceToVolumeRegistration>::Create(
        new EbbRTSliceToVolumeRegistration(_reconstructor, kPrinterEbbId),
        kPrinterEbbId);

    // returns a future for the EbbRef in AppMain.cc
    return ebbrt::global_id_map
        ->Set(kPrinterEbbId, ebbrt::messenger->LocalNetworkId().ToBytes())
        // the Future<void> f is returned by the ebbrt::global_id_map
        .Then([ebbref](ebbrt::Future<void> f) {
          // f.Get() ensures the gobal_id_map has completed
          f.Get();
          return ebbref;
        });
  }

  // this uses the void promise to return a future that will
  // get invoked when mypromose gets set with some value
  ebbrt::Future<void> waitReceive() { return mypromise.GetFuture(); }

  void Print(ebbrt::Messenger::NetworkId nid, const char* str);

  void ReceiveMessage(ebbrt::Messenger::NetworkId nid,
                      std::unique_ptr<ebbrt::IOBuf>&& buffer);

  EbbId getEbbId() { return ebbid; }

  void initMatrix(irtkMatrix& m, double val, int row, int col) {
    for (int i = 0; i < row; i++) {
      for (int j = 0; j < col; j++) {
        m.Put(i, j, val);
      }
    }
  }

  void runJob();
  void slicetovolume();
  
  void destroy() { delete this; }
};

#endif
