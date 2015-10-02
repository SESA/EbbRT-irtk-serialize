//          Copyright Boston University SESA Group 2013 - 2014.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <signal.h>

#include <boost/filesystem.hpp>

#include <ebbrt/Context.h>
#include <ebbrt/ContextActivation.h>
#include <ebbrt/GlobalIdMap.h>
#include <ebbrt/StaticIds.h>
#include <ebbrt/NodeAllocator.h>
#include <ebbrt/Runtime.h>

#include "Printer.h"
#include "EbbRTStackRegistrations.h"

int main(int argc, char** argv) {
  auto bindir = boost::filesystem::system_complete(argv[0]).parent_path() /
                "/bm/AppMain.elf32";

  ebbrt::Runtime runtime;
  ebbrt::Context c(runtime);
  ebbrt::ContextActivation activation(c);

  irtkReconstruction* reconstructor;
  vector<irtkRealImage> stacks;
  vector<irtkRigidTransformation> stack_transformation;
  int templateNumber = 0;
  irtkGreyImage target;
  irtkRigidTransformation offset;
  bool useExternalTarget = false;

  /*
   * event_manager -> Spawn(..) puts the code block on the event queue
   * to be executed after resolving the Future
   *
   */
  ebbrt::event_manager
      ->Spawn([&reconstructor, &stacks, &stack_transformation, templateNumber,
               &target, &offset, useExternalTarget, bindir]() {

        // EbbRTStackRegistrations::Create(..) returns a
        // Future<EbbRTStackRegistrationsEbbRef>
        // that gets accessed on the Then(..) call, f.Get() ensures the EbbRef
        // was created successfully, then we can allocate the baremetal node on
        // the backend.
        EbbRTStackRegistrations::Create(reconstructor, stacks,
                                        stack_transformation, templateNumber,
                                        target, offset, useExternalTarget)
            // Then(...) gets the EbbRef
            .Then([bindir](ebbrt::Future<EbbRTStackRegistrationsEbbRef> f) {
              // ensures it was created
              EbbRTStackRegistrationsEbbRef ref = f.Get();

              // allocated baremetal AppMain.elf32
              ebbrt::node_allocator->AllocateNode(bindir.string());

              // test code to get EbbId
              std::cout << "EbbId: " << ref->getEbbId() << std::endl;
            });
      });

  c.Run();

  return 0;
}
