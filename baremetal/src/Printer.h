//          Copyright Boston University SESA Group 2013 - 2014.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
#ifndef APPS_PINGPONG_BAREMETAL_SRC_PRINTER_H_
#define APPS_PINGPONG_BAREMETAL_SRC_PRINTER_H_

#include <string>
#include <string.h>

#include <ebbrt/Message.h>

#include <irtkReconstructionGPU.h>
#include <irtkResampling.h>
#include <irtkRegistration.h>
#include <irtkImageRigidRegistration.h>
#include <irtkImageRigidRegistrationWithPadding.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>

#include "../../src/StaticEbbIds.h"

class Printer : public ebbrt::Messagable<Printer> {
 public:
  explicit Printer(ebbrt::Messenger::NetworkId nid);

  static Printer& HandleFault(ebbrt::EbbId id);

  void Print(const char* string);
  void ReceiveMessage(ebbrt::Messenger::NetworkId nid,
                      std::unique_ptr<ebbrt::IOBuf>&& buffer);

 private:
  ebbrt::Messenger::NetworkId remote_nid_;
};

constexpr auto printer = ebbrt::EbbRef<Printer>(kPrinterEbbId);

/*struct POINT3D
{
  short x;
  short y;
  short z;
  float value;
};

typedef std::vector<POINT3D> VOXELCOEFFS;
typedef std::vector<std::vector<VOXELCOEFFS> > SLICECOEFFS;
*/


#endif  // APPS_HELLOWORLD_BAREMETAL_SRC_PRINTER_H_
