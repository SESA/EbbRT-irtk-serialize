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

//EBBRT_PUBLISH_TYPE(, EbbRTStackRegistrations);

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
  std::string bmnid;
  irtkMatrix mym;
  // this is used to save and load context
  ebbrt::EventManager::EventContext* emec{nullptr};

  // a void promise where we want to start some computation after some other
  // function has finished
  ebbrt::Promise<void> mypromise;

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
    bmnid = "";
  }

  //Default Constructor
  //EbbRTStackRegistrations() : Messagable<EbbRTStackRegistrations>(0);
  
  /* static ebbrt::Future<EbbRTStackRegistrationsEbbRef> */
  /* Create(irtkReconstruction* _reconstructor) { */
  /*   // calls the constructor in this class and returns a EbbRef */
  /*   auto ebbref = SharedEbb<EbbRTStackRegistrations>::Create( */
  /*       new EbbRTStackRegistrations(_reconstructor), kPrinterEbbId); */

  /*   // returns a future for the EbbRef in AppMain.cc */
  /*   return ebbrt::global_id_map */
  /*       ->Set(kPrinterEbbId, ebbrt::messenger->LocalNetworkId().ToBytes()) */
  /*       // the Future<void> f is returned by the ebbrt::global_id_map */
  /*       .Then([ebbref](ebbrt::Future<void> f) { */
  /*         // f.Get() ensures the gobal_id_map has completed */
  /*         f.Get(); */
  /*         return ebbref; */
  /*       }); */
  /* } */

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

  // this uses the void promise to return a future that will
  // get invoked when mypromose gets set with some value
  ebbrt::Future<void> waitReceive() { return mypromise.GetFuture(); }

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
    std::cout << "Received msg: " << nid.ToString() << ": " << output << "\n";
    // std::this_thread::sleep_for(std::chrono::seconds(1));

    // serialize a vector
    /*std::vector<int> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa & v;*/
    /*
    if (output == "ping") {
      std::cout << nid.ToString() << ": " << output << "\n";
      Print(nid, ofs.str().c_str());
    } else if (strcmp(output.c_str(), "pong") == 0) {
      std::cout << nid.ToString() << ": " << output << "\n";
      Print(nid, "ping");
    }
    */

    // stores the received nid as a string of bytes
    if (bmnid.empty()) {
      bmnid = nid.ToBytes();
      std::cout << "isEmpty()? " << bmnid.empty() << std::endl;

      //*****
      // this SetValue() sets up the promise so that the waitReceive()
      // future is fulfilled and starts the next computation in AppMain()
      // this call ensures that we don't call SendMessage until we know
      // the bm node has been initialized
      mypromise.SetValue();
    }

    if (output[0] == 'A') {
      std::cout << "Parsing it back" << std::endl;

      std::istringstream ifs;
      ifs.str(output.substr(2, output.length() - 2));
      boost::archive::text_iarchive ia(ifs);

      ia& mym;

      // reactivate context in runJob2()
      ebbrt::event_manager->ActivateContext(std::move(*emec));

      // std::cout << "Received matrix: " << std::endl;
      // mym.Print();
      // ebbrt::active_context->io_service_.stop();
    } else if (output[0] == 'B') {
      std::cout << "Parsing it back" << std::endl;

      std::istringstream ifs;
      ifs.str(output.substr(2, output.length() - 2));
      boost::archive::text_iarchive ia(ifs);

      ia& mym;

      // reactivate context in runJob2()
      ebbrt::event_manager->ActivateContext(std::move(*emec));
    }

    // ebbrt::Messenger::NetworkId bmnidd;
    // try to reintepret it back as a nid
    // Print(bmnidd.FromBytes(reinterpret_cast<const unsigned
    // char*>(bmnid.c_str()), bmnid.size()), ofs.str().c_str());

    // wait for 1 second before shutting down, this doesn't use Ctrl-C
    // std::this_thread::sleep_for(std::chrono::seconds(1));
    // std::cout << "Shutting down ..." << std::endl;

    // referencing active_context defined in
    // EbbRT/hosted/src/include/ebbrt/Context.h
    // gets this current ebbs active context
    // ebbrt::active_context->io_service_.stop();
  }

  EbbId getEbbId() { return ebbid; }

  void initMatrix(irtkMatrix& m, double val, int row, int col) {
    for (int i = 0; i < row; i++) {
      for (int j = 0; j < col; j++) {
        m.Put(i, j, val);
      }
    }
  }

  void slicetovolume() {}

  void runJob3() {
    // get the event manager context and save it
    ebbrt::EventManager::EventContext context;

    irtkMatrix m1(5, 5);

    initMatrix(m1, 1.0, 5, 5);

    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa& m1;
    std::string ts = "B " + ofs.str();

    // sending serialized vector to bm
    ebbrt::Messenger::NetworkId bmnidd;
    Print(
        bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
                         bmnid.size()),
        ts.c_str());
    std::cout << "Sent vector " << std::endl;

    emec = &context;
    ebbrt::event_manager->SaveContext(*emec);

    // this should execute when the saved context is loaded again
    std::cout << "Received matrix: " << std::endl;
    mym.Print();
    ebbrt::active_context->io_service_.stop();
  }

  void runJob2() {
    // get the event manager context and save it
    ebbrt::EventManager::EventContext context;

    irtkMatrix m1(2, 2);
    irtkMatrix m2(2, 2);
    // irtkMatrix m3;

    initMatrix(m1, 2.0, 2, 2);
    initMatrix(m2, 3.0, 2, 2);

    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa& m1& m2;

    std::string ts = "A " + ofs.str();
    // std::cout << ofs.str() << std::endl;
    // std::cout << ts << std::endl;

    // m3 = m1 * m2;
    // m3.Print();

    // sending serialized vector to bm
    ebbrt::Messenger::NetworkId bmnidd;
    Print(
        bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
                         bmnid.size()),
        ts.c_str());
    std::cout << "Sent vector " << std::endl;

    emec = &context;
    ebbrt::event_manager->SaveContext(*emec);

    // this should execute when the saved context is loaded again
    std::cout << "Received matrix: " << std::endl;
    mym.Print();
    ebbrt::active_context->io_service_.stop();
  }

  void runJob(int size) {
    std::cout << " in Ebbrt" << std::endl;
    // get the event manager context and save it
    ebbrt::EventManager::EventContext context;

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

      // serialize mo and m
      ebbrt::Messenger::NetworkId bmnidd;
      std::ostringstream ofs;
      boost::archive::text_oarchive oa(ofs);
      oa& mo& m;

      std::string ts = "A " + ofs.str();
      Print(bmnidd.FromBytes(
                reinterpret_cast<const unsigned char*>(bmnid.c_str()),
                bmnid.size()),
            ts.c_str());
      std::cout << "Sent vector " << std::endl;

      emec = &context;
      ebbrt::event_manager->SaveContext(*emec);

      m = mym;
      // m = m * mo;
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
