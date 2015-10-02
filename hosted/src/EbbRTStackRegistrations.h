#ifndef _EBBRT_STACK_REGISTRATIONS_H_
#define _EBBRT_STACK_REGISTRATIONS_H_

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

class EbbRTStackRegistrations;
typedef EbbRef<EbbRTStackRegistrations> EbbRTStackRegistrationsEbbRef;

class EbbRTStackRegistrations : public SharedEbb<EbbRTStackRegistrations>,
                                public Messagable<EbbRTStackRegistrations> {

  irtkReconstruction* reconstructor;
  vector<irtkRealImage>& stacks;
  vector<irtkRigidTransformation>& stack_transformations;
  int templateNumber;
  irtkGreyImage& target;
  irtkRigidTransformation& offset;
  bool externalTemplate;
  EbbId ebbid;

 public:
  EbbRTStackRegistrations(
      irtkReconstruction* _reconstructor, vector<irtkRealImage>& _stacks,
      vector<irtkRigidTransformation>& _stack_transformations,
      int _templateNumber, irtkGreyImage& _target,
      irtkRigidTransformation& _offset, EbbId _ebbid,
      bool _externalTemplate = false)
      : Messagable<EbbRTStackRegistrations>(_ebbid),
        reconstructor(_reconstructor), stacks(_stacks),
        stack_transformations(_stack_transformations), target(_target),
        offset(_offset) {
    templateNumber = _templateNumber;
    externalTemplate = _externalTemplate;
    ebbid = _ebbid;
  }

  // Create function that returns a future ebbref
  static ebbrt::Future<EbbRTStackRegistrationsEbbRef>
  Create(irtkReconstruction* _reconstructor, vector<irtkRealImage>& _stacks,
         vector<irtkRigidTransformation>& _stack_transformations,
         int _templateNumber, irtkGreyImage& _target,
         irtkRigidTransformation& _offset, bool _externalTemplate = false) {
    // was taken out since not sure how to pass the locally allocated ebbid to
    // the backend,
    // perhaps through the global id map?
    //   EbbId id = ebb_allocator->AllocateLocal()) {

    // calls the constructor in this class and returns a EbbRef
    auto ebbref = SharedEbb<EbbRTStackRegistrations>::Create(
        new EbbRTStackRegistrations(
            _reconstructor, _stacks, _stack_transformations, _templateNumber,
            _target, _offset, kPrinterEbbId, _externalTemplate),
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

  void Print(ebbrt::Messenger::NetworkId nid, const char* str) {
    auto len = strlen(str) + 1;
    auto buf = ebbrt::MakeUniqueIOBuf(len);
    snprintf(reinterpret_cast<char*>(buf->MutData()), len, "%s", str);
    SendMessage(nid, std::move(buf));
  }

  /*void ReceiveMessage(Messenger::NetworkId nid,
                      unique_ptr<ebbrt::IOBuf>&& buffer) {
    return;
    }*/

  void ReceiveMessage(ebbrt::Messenger::NetworkId nid,
                      std::unique_ptr<ebbrt::IOBuf>&& buffer) {
    /***********************************
     * std::string(reinterpret_cast<const char*> buffer->Data(),
     *buffer->Length())
     *
     * Using the above code, I had to first do strcmp(output.c_str(), ...) to
     *ensure
     * it matched the input string.
     * Direct comparison using "==" seems to be working when I don't pass in the
     *Length() as second arg
     **********************************/
    auto output = std::string(reinterpret_cast<const char*>(buffer->Data()));
    std::this_thread::sleep_for(std::chrono::seconds(1));

    // serialize a vector
    std::vector<int> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa& v;

    if (output == "ping") {
      std::cout << nid.ToString() << ": " << output << "\n";
      Print(nid, ofs.str().c_str());
    } else if (strcmp(output.c_str(), "pong") == 0) {
      std::cout << nid.ToString() << ": " << output << "\n";
      Print(nid, "ping");
    }

    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "Shutting down ..." << std::endl;

    // referencing active_context defined in
    // EbbRT/hosted/src/include/ebbrt/Context.h
    // gets this current ebbs active context
    ebbrt::active_context->io_service_.stop();
  }

  EbbId getEbbId() { return ebbid; }

  void runJob(int size) {
    std::cout << " in Ebbrt" << std::endl;

    for (int i = 0; i < size; ++i) {
      // do not perform registration for template
      if (i == templateNumber)
        continue;

      // rigid registration object
      irtkImageRigidRegistrationWithPadding registration;

      // set target and source (need to be converted to irtkGreyImage)
      irtkGreyImage source = stacks[i];

      // include offset in trasformation
      irtkMatrix mo = offset.GetMatrix();
      irtkMatrix m = stack_transformations[i].GetMatrix();
      m = m * mo;
      stack_transformations[i].PutMatrix(m);

      // perform rigid registration
      registration.SetInput(&target, &source);
      registration.SetOutput(&stack_transformations[i]);
      if (externalTemplate) {
        registration.GuessParameterThickSlicesNMI();
      } else {
        registration.GuessParameterThickSlices();
      }
      registration.SetTargetPadding(0);
      registration.Run();

      mo.Invert();
      m = stack_transformations[i].GetMatrix();
      m = m * mo;
      stack_transformations[i].PutMatrix(m);
    }
  }

  void operator()(const blocked_range<size_t>& r) const;
  // execute
  void operator()() const;

  void destroy() { delete this; }
};

#endif
