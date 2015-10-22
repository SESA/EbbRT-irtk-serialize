#include "EbbRTCoeffInit.h"

EBBRT_PUBLISH_TYPE(, EbbRTCoeffInit);

void EbbRTCoeffInit::Print(ebbrt::Messenger::NetworkId nid,
                                           const char* str) {
  auto len = strlen(str) + 1;
  auto buf = ebbrt::MakeUniqueIOBuf(len);
  snprintf(reinterpret_cast<char*>(buf->MutData()), len, "%s", str);
  SendMessage(nid, std::move(buf));
}

void EbbRTCoeffInit::ReceiveMessage(
    ebbrt::Messenger::NetworkId nid, std::unique_ptr<ebbrt::IOBuf>&& buffer) {
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
  } 
  else if (output[0] == 'C') 
  {
      std::cout << "Parsing it back" << std::endl;

      std::istringstream ifs;
      ifs.str(output.substr(2, output.length() - 2));
      boost::archive::text_iarchive ia(ifs);

      ia & _slice;// & _transformations & slicecoeffs & slice_inside;
            
      
      // reactivate context in runJob2()
      ebbrt::event_manager->ActivateContext(std::move(*emec));
    }
}

void EbbRTCoeffInit::runJob() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;
  int i = 0;
  
  for(i=0;i<2;i++)
  {
      irtkMatrix m1(2, 2);
      irtkMatrix m2(2, 2);
      // irtkMatrix m3;

      initMatrix(m1, 2.0, 2, 2);
      initMatrix(m2, 3.0, 2, 2);
      irtkRealImage test;
      
      std::cout << "before" <<std::endl;
      
      std::ostringstream ofs;
      boost::archive::text_oarchive oa(ofs);
      oa & m1 & m2 & test;
      
      std::cout << "after" <<std::endl;

      std::string ts = "A " + ofs.str();
      
      // sending serialized vector to bm
      ebbrt::Messenger::NetworkId bmnidd;
      Print(bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
			     bmnid.size()),
	    ts.c_str());
      std::cout << "Sent vector " << std::endl;
      
      emec = &context;
      ebbrt::event_manager->SaveContext(*emec);
      
      // this should execute when the saved context is loaded again
      std::cout << "Received matrix: " << std::endl;
      mym.Print();
  }

  ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::coeffinit() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;

  irtkImageAttributes attr = reconstructor->_reconstructed.GetImageAttributes();
  size_t max_slices = reconstructor->_slices.size();
  size_t inputIndex = 0;

  std::cout << "Inside coeffinit() " << max_slices 
	    << " " << inputIndex << std::endl;

  for (inputIndex = 0; inputIndex != max_slices; inputIndex++) 
  {
      // serialization
      ebbrt::Messenger::NetworkId bmnidd;
      std::ostringstream ofs;
      boost::archive::text_oarchive oa(ofs);

      std::cout << "Before serialize" << std::endl;
      
      oa & reconstructor->_slices[inputIndex]
	  & reconstructor->_reconstructed
	  & reconstructor->_quality_factor 
	  & reconstructor->_transformations[inputIndex]
	  & reconstructor->_mask;
      
      std::cout << "after serialize" << std::endl;
      
      std::string ts = "C " + ofs.str();
      Print(
	  bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
			   bmnid.size()),
	  ts.c_str());

      std::cout << "Sending data" << std::endl;
      emec = &context;
      ebbrt::event_manager->SaveContext(*emec);

      std::cout << "Received back " <<std::endl;
      reconstructor->_slices[inputIndex] = _slice;
      //reconstructor->_transformations[inputIndex] = _transformations;
      //reconstructor->_volcoeffs[inputIndex] = slicecoeffs;
      //reconstructor->_slice_inside_cpu[inputIndex] = slice_inside;
      
      std::cout << "EbbRTCoeffInit done " << std::endl;
  }

  ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::coeffinit2() {

    std::cout << "Inside coeffinit" << std::endl;
    size_t inputIndex = 0;
    
    for (inputIndex=0;inputIndex!=reconstructor->_slices.size();++inputIndex) 
    {
	bool slice_inside;
	
	// get resolution of the volume
	double vx, vy, vz;
	reconstructor->_reconstructed.GetPixelSize(&vx, &vy, &vz);
	// volume is always isotropic
	double res = vx;
	
	// start of a loop for a slice inputIndex
	cout << inputIndex << " ";
	
	// read the slice
	irtkRealImage &slice = reconstructor->_slices[inputIndex];
	
	// prepare structures for storage
	POINT3D p;
	VOXELCOEFFS empty;
	SLICECOEFFS slicecoeffs(slice.GetX(),
				vector<VOXELCOEFFS>(slice.GetY(), empty));

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
	double size = res / reconstructor->_quality_factor;

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
		    reconstructor->_transformations[inputIndex].Transform(x, y, z);
		    reconstructor->_reconstructed.WorldToImage(x, y, z);
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
				reconstructor->_transformations[inputIndex].Transform(x, y,
										      z);
				// Change to image coordinates
				reconstructor->_reconstructed.WorldToImage(x, y, z);

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
				    if ((l >= 0) && (l < reconstructor->_reconstructed.GetX()))
					for (m = ny; m <= ny + 1; m++)
					    if ((m >= 0) &&
						(m < reconstructor->_reconstructed.GetY()))
						for (n = nz; n <= nz + 1; n++)
						    if ((n >= 0) &&
							(n < reconstructor->_reconstructed.GetZ())) {
							weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) *
							    (1 - fabs(n - z));
							sum += weight;
							if (reconstructor->_mask(l, m, n) == 1) {
							    inside = true;
							    slice_inside = true;
							}
						    }
				// if there were no voxels do nothing
				if ((sum <= 0) || (!inside))
				    continue;
				// now calculate the transformed PSF
				for (l = nx; l <= nx + 1; l++)
				    if ((l >= 0) && (l < reconstructor->_reconstructed.GetX()))
					for (m = ny; m <= ny + 1; m++)
					    if ((m >= 0) &&
						(m < reconstructor->_reconstructed.GetY()))
						for (n = nz; n <= nz + 1; n++)
						    if ((n >= 0) &&
							(n < reconstructor->_reconstructed.GetZ())) {
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
      
	reconstructor->_volcoeffs[inputIndex] = slicecoeffs;
	reconstructor->_slice_inside_cpu[inputIndex] = slice_inside;

    } // end of loop through the slices

    ebbrt::active_context->io_service_.stop();
}
