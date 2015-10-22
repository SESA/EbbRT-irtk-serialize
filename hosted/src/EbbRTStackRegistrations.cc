#include "EbbRTStackRegistrations.h"

EBBRT_PUBLISH_TYPE(, EbbRTStackRegistrations);


/********EbbRTStackRegistrations::EbbRTStackRegistrations()
{
    reconstructor = NULL;
    stacks;
    stack_transformations = NULL;
    offset = NULL;
    target = NULL;
    templateNUmber = 0;
    ebbid = 0;
    externalTemplate = false;
    bmnid = "";
    //Messagable<EbbRTStackRegistrations>(0);
    }*/


// void EbbRTStackRegistrations::operator()(const blocked_range<size_t>& r) const {
//   for (int i = r.begin(); i != r.end(); ++i) {

//     // do not perform registration for template
//     if (i == templateNumber)
//       continue;

//     // test
//     // rigid registration object
//     irtkImageRigidRegistrationWithPadding registration;
//     // irtkRigidTransformation transformation = stack_transformations[i];

//     // set target and source (need to be converted to irtkGreyImage)
//     irtkGreyImage source = stacks[i];

//     // include offset in trasformation
//     irtkMatrix mo = offset.GetMatrix();
//     irtkMatrix m = stack_transformations[i].GetMatrix();
//     m = m * mo;
//     stack_transformations[i].PutMatrix(m);

//     // perform rigid registration
//     registration.SetInput(&target, &source);
//     registration.SetOutput(&stack_transformations[i]);
//     if (externalTemplate) {
//       registration.GuessParameterThickSlicesNMI();
//     } else {
//       registration.GuessParameterThickSlices();
//     }
//     registration.SetTargetPadding(0);
//     registration.Run();

//     mo.Invert();
//     m = stack_transformations[i].GetMatrix();
//     m = m * mo;
//     stack_transformations[i].PutMatrix(m);

//     // stack_transformations[i] = transformation;

//     // save volumetric registrations
//     /*if (reconstructor->_debug) {
//       //buffer to create the name
//       char buffer[256];
//       registration.irtkImageRegistration::Write((char *) "parout-volume.rreg");
//       sprintf(buffer, "stack-transformation%i.dof.gz", i);
//       stack_transformations[i].irtkTransformation::Write(buffer);
//       target.Write("target.nii.gz");
//       sprintf(buffer, "stack%i.nii.gz", i);
//       stacks[i].Write(buffer);
//       }*/
//   }
// }

// // execute
// void EbbRTStackRegistrations::operator()() const {
//   // task_scheduler_init init(tbb_no_threads);
//   //task_scheduler_init init(1);
//   //parallel_for(blocked_range<size_t>(0, stacks.size()), *this);
//   //init.terminate();
// }
