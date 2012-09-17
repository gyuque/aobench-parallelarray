/*   ooo ooo ooo ooo ooo   aobench   ooo ooo ooo ooo ooo
 *
 *            Original version by Syoyo Fujita
 *            http://code.google.com/p/aobench/
 *
 * Javascript + River Trail version by Satoshi Ueyama, 2012
 */

(function() {
	//'use strict';
	// Unfortunately, River Trail fails in strict mode

	// Function shortcuts
	var sqrt = Math.sqrt;
	var fabs = Math.abs;
	var drand48 = Math.random;
	var M_PI = Math.PI;
	var cos = Math.cos;
	var sin = Math.sin;
	
	var WIDTH  = 128;
	var HEIGHT = 128;
	var NSUBSAMPLES = 2;
	var NAO_SAMPLES = 8;
	
	var gScene = null;
	
	// Geometries ---------------------------------
	function Vec() {
		this.x = 0;
		this.y = 0;
		this.z = 0;
	}
	
	function vdot(v0, v1) {
		return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
	}
	
	Vec.prototype = {
		normalize: function() {
			var length = sqrt(vdot(this, this));

			if (fabs(length) > 1.0e-17) {
				this.x /= length;
				this.y /= length;
				this.z /= length;
			}
		},
		
		cross: function(v0, v1) {
			this.x = v0.y * v1.z - v0.z * v1.y;
			this.y = v0.z * v1.x - v0.x * v1.z;
			this.z = v0.x * v1.y - v0.y * v1.x;
		}
	};
	
	function Sphere() {
		this.center = new Vec();
		this.radius = 0.0;
	}
	
	function Plane() {
		this.p = new Vec(); // position
		this.n = new Vec(); // normal
	}
	
	// Rendering informations -----------------------
	function Ray() {
		this.org = new Vec();
		this.dir = new Vec();
	}
	
	function Isect() {
		this.t = 0; //distance
		this.p = new Vec(); // position
		this.n = new Vec(); // normal
		this.hit = false; // hit flag
	};
	
	// Hit tests
	function raySphereIntersect(isect, ray, sphere, tempVec) {
		var rs = tempVec;
		
		rs.x = ray.org.x - sphere.center.x;
		rs.y = ray.org.y - sphere.center.y;
		rs.z = ray.org.z - sphere.center.z;
		var B = vdot(rs, ray.dir);
		var C = vdot(rs, rs) - sphere.radius * sphere.radius;
		var D = B * B - C;

		if (D > 0.0) {
			var t = -B - sqrt(D);
			
			if ((t > 0.0) && (t < isect.t)) {
				isect.t = t;
				isect.hit = true;

				isect.p.x = ray.org.x + ray.dir.x * t;
				isect.p.y = ray.org.y + ray.dir.y * t;
				isect.p.z = ray.org.z + ray.dir.z * t;

				isect.n.x = isect.p.x - sphere.center.x;
				isect.n.y = isect.p.y - sphere.center.y;
				isect.n.z = isect.p.z - sphere.center.z;

				isect.n.normalize();
			}
		}
	}
	
	function rayPlaneIntersect(isect, ray, plane)
	{
		var d = -vdot(plane.p, plane.n);
		var v = vdot(ray.dir, plane.n);

		if (fabs(v) < 1.0e-17) return;

		var t = -(vdot(ray.org, plane.n) + d) / v;

		if ((t > 0.0) && (t < isect.t)) {
			isect.t = t;
			isect.hit = true;

			isect.p.x = ray.org.x + ray.dir.x * t;
			isect.p.y = ray.org.y + ray.dir.y * t;
			isect.p.z = ray.org.z + ray.dir.z * t;

			isect.n.x = plane.n.x;
			isect.n.y = plane.n.y;
			isect.n.z = plane.n.z;
		}
	}

	// Renderer functions

	function orthoBasis(basis, n)
	{
		basis[2] = n;
		basis[1].x = 0.0; basis[1].y = 0.0; basis[1].z = 0.0;

		if ((n.x < 0.6) && (n.x > -0.6)) {
			basis[1].x = 1.0;
	    } else if ((n.y < 0.6) && (n.y > -0.6)) {
			basis[1].y = 1.0;
	    } else if ((n.z < 0.6) && (n.z > -0.6)) {
			basis[1].z = 1.0;
	    } else {
			basis[1].x = 1.0;
	    }

	    basis[0].cross(basis[1], basis[2]);
	    basis[0].normalize();

	    basis[1].cross(basis[2], basis[0]);
	    basis[1].normalize();
	}


	function ambientOcclusion(scn, col, isect, rsTempVec, aoRayTemp, aoOrgTemp, aoIsectTemp, basisTemp) {
		var i, j;
		var ntheta = NAO_SAMPLES;
		var nphi   = NAO_SAMPLES;
		var eps = 0.0001;

		var p = aoOrgTemp;
		p.x = isect.p.x + eps * isect.n.x;
		p.y = isect.p.y + eps * isect.n.y;
		p.z = isect.p.z + eps * isect.n.z;
		
		var basis = basisTemp;
		orthoBasis(basis, isect.n);
		
		var occlusion = 0.0;
		var pi2 = M_PI * 2.0;

		// Do monte carlo sampling for secondary rays
		for (j = 0; j < ntheta; ++j) {
			for (i = 0; i < nphi; ++i) {
				var theta = sqrt(drand48());
				var phi   = pi2 * drand48();

				// Select a random direction
				var x = cos(phi) * theta;
				var y = sin(phi) * theta;
				var z = sqrt(1.0 - theta * theta);
				
				// Transform ray direction on local plane to global coordinate
				var rx = x * basis[0].x + y * basis[1].x + z * basis[2].x;
				var ry = x * basis[0].y + y * basis[1].y + z * basis[2].y;
				var rz = x * basis[0].z + y * basis[1].z + z * basis[2].z;
				
				var ray = aoRayTemp;
				ray.org.x = p.x;
				ray.org.y = p.y;
				ray.org.z = p.z;
				ray.dir.x = rx;
				ray.dir.y = ry;
				ray.dir.z = rz;
				
				var occIsect = aoIsectTemp;
				occIsect.t   = 1.0e+17;
				occIsect.hit = false;

				// Cast a ray
				raySphereIntersect(occIsect, ray, scn.spheres[0], rsTempVec);
				raySphereIntersect(occIsect, ray, scn.spheres[1], rsTempVec);
				raySphereIntersect(occIsect, ray, scn.spheres[2], rsTempVec);
				rayPlaneIntersect (occIsect, ray, scn.plane);
				if (occIsect.hit) occlusion += 1.0;
			}
		}
		
		occlusion = (ntheta * nphi - occlusion) / (ntheta * nphi);

		// Final irradiance
		col.x = occlusion;
		col.y = occlusion;
		col.z = occlusion;
	}	
	
	function initScene()
	{
		var scn = {
			spheres: [new Sphere(), new Sphere(), new Sphere()],
			plane: new Plane()
		};
		
		var spheres = scn.spheres;
		var plane = scn.plane;
		
	    spheres[0].center.x = -2.0;
	    spheres[0].center.y =  0.0;
	    spheres[0].center.z = -3.5;
	    spheres[0].radius = 0.5;

	    spheres[1].center.x = -0.5;
	    spheres[1].center.y =  0.0;
	    spheres[1].center.z = -3.0;
	    spheres[1].radius = 0.5;

	    spheres[2].center.x =  1.0;
	    spheres[2].center.y =  0.0;
	    spheres[2].center.z = -2.2;
	    spheres[2].radius = 0.5;

	    plane.p.x = 0.0;
	    plane.p.y = -0.5;
	    plane.p.z = 0.0;

	    plane.n.x = 0.0;
	    plane.n.y = 1.0;
	    plane.n.z = 0.0;
	
		return scn;
	}
	
	function render(img, w, h, nsubsamples)
	{
		// Work areas (to suppress mass allocation)
		var rsTempVec = new Vec();
		var aoRayTemp = new Ray();
		var aoOrgTemp = new Vec();
		var aoIsectTemp = new Isect();
		var basisTemp = [new Vec(), new Vec(), new Vec()];

		var col = new Vec();
	    var x, y;
	    var u, v;

		var ray = new Ray();
		var isect = new Isect();
	    var fimg = allocateFloats(w * h * 3);

		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				for (v = 0; v < nsubsamples; ++v) {
	                for (u = 0; u < nsubsamples; ++u) {
						// Make normalized coordinate
						var px = (x + (u / nsubsamples) - (w / 2.0)) / (w / 2.0);
						var py = -(y + (v / nsubsamples) - (h / 2.0)) / (h / 2.0);
						
						// Eye position
						ray.org.x = 0.0;
						ray.org.y = 0.0;
						ray.org.z = 0.0;
						
						// Ray direction from the eye
						ray.dir.x = px;
						ray.dir.y = py;
						ray.dir.z = -1.0;
						ray.dir.normalize();
						
						isect.t   = 1.0e+17;
						isect.hit = false;
						
						// Cast the primary ray
						raySphereIntersect(isect, ray, gScene.spheres[0], rsTempVec);
						raySphereIntersect(isect, ray, gScene.spheres[1], rsTempVec);
						raySphereIntersect(isect, ray, gScene.spheres[2], rsTempVec);
						rayPlaneIntersect (isect, ray, gScene.plane);

						if (isect.hit) {
							// Do secondary ray tracing
							ambientOcclusion(gScene, col, isect, rsTempVec, aoRayTemp, aoOrgTemp, aoIsectTemp, basisTemp);
							fimg[3 * (y * w + x) + 0] += col.x;
							fimg[3 * (y * w + x) + 1] += col.y;
							fimg[3 * (y * w + x) + 2] += col.z;
						}
					}
				}
				
				var bufpos = 3 * (y * w + x);
				fimg[bufpos + 0] /= nsubsamples * nsubsamples;
				fimg[bufpos + 1] /= nsubsamples * nsubsamples;
				fimg[bufpos + 2] /= nsubsamples * nsubsamples;
				
				img[bufpos + 0] = fimg[bufpos + 0] * 255;
				img[bufpos + 1] = fimg[bufpos + 1] * 255;
				img[bufpos + 2] = fimg[bufpos + 2] * 255;
			}
		}
	}
	
	function emitToCanvas(g, img, w, h) {
		var x, y;
		var imageData = g.getImageData(0, 0, w, h);
		var pixs = imageData.data;
		
		var rpos = 0;
		var wpos = 0;
		for (y = 0;y < h;++y) {
			for (x = 0;x < w;++x) {
				pixs[wpos++] = img[rpos++];
				pixs[wpos++] = img[rpos++];
				pixs[wpos++] = img[rpos++];
				pixs[wpos++] = 255;
			}
		}
		
		g.putImageData(imageData, 0, 0);
	}
	
	// Main
	window.launch = function() {
		var outCanvas = document.getElementById("out-canvas");
		outCanvas.width  = WIDTH;
		outCanvas.height = HEIGHT;
		/*
		var testpa = new ParallelArray([0,1,2,3,4,5,6,7]);
		testpa.shape = [2, 4];
		testpa.strides = [4, 1];
		console.log( testpa.combine(function(i) {
			var topIndex = i[0] - 0;
			var a = this.get(topIndex);
			return a[0] + a[2] * Math.sqrt(Math.cos(1));
		}) );
		
		console.log(new ParallelArray([ [0,1,2,3],[3,4,5,6] ]));
		*/
		
		function startAO(par) {
			var g = outCanvas.getContext("2d");
		
			var img = allocateBytes(WIDTH * HEIGHT * 3);
			gScene = initScene();
			
			if (!par)
				render(img, WIDTH, HEIGHT, NSUBSAMPLES);
			else
				parRender(img, WIDTH, HEIGHT, NSUBSAMPLES);
		
			emitToCanvas(g, img, WIDTH, HEIGHT);
		}
		
		var btn = document.getElementById("start-btn");
		var par_btn = document.getElementById("par-start-btn");
		
		function buttonHandler(par) {
			setTimeout(function() {
				var startT = new Date();
				startAO(par);
				
				var rout = document.getElementById("result-out");
				rout.innerHTML = "Parallel: " + (par ? "Yes" : "No") + "<br>" +(new Date() - startT) + " ms - " + 
				                 WIDTH+"x"+HEIGHT+" px, "+NSUBSAMPLES+"x"+NSUBSAMPLES+" subsamples";
				
				setTimeout(function() {
					rout.style.opacity = 1;
					rout.style.top = "10px";
				}, 100);
			}, 50);

			btn.style.display = "none";
			par_btn.style.display = "none";
		}
		
		btn.addEventListener('click', function() {
			buttonHandler(false);
		}, false);

		par_btn.addEventListener('click', function() {
			buttonHandler(true);
		}, false);
	};
	
	// Utilities
	function allocateBytes(size) {
		var a = new Uint8ClampedArray(size);
		return a;
	}

	function allocateFloats(size) {
		var a = new Float32Array(size);
		return a;
	}
	
	// ++++ Parallel version ++++
	function addAOEntry(outArray, ai, scn, x, y, isect, aoOrgTemp, basisTemp) {
		var i, j;
		var ntheta = NAO_SAMPLES;
		var nphi   = NAO_SAMPLES;
		var eps = 0.0001;

		var p = aoOrgTemp;
		p.x = isect.p.x + eps * isect.n.x;
		p.y = isect.p.y + eps * isect.n.y;
		p.z = isect.p.z + eps * isect.n.z;
		
		var basis = basisTemp;
		orthoBasis(basis, isect.n);
		
		var occlusion = 0.0;
		var pi2 = M_PI * 2.0;

		var sphereList = scn.spheres;
		// Do monte carlo sampling for secondary rays
		for (j = 0; j < ntheta; ++j) {
			for (i = 0; i < nphi; ++i) {
				var theta = sqrt(drand48());
				var phi   = pi2 * drand48();
				for (var si = 0;si < 3;++si) {
					var sphere = sphereList[si];
					
					outArray[ai++] = x;
					outArray[ai++] = y;
					outArray[ai++] = theta;
					outArray[ai++] = phi;

					outArray[ai++] = basis[1].x;
					outArray[ai++] = basis[1].y;
					outArray[ai++] = basis[1].z;

					outArray[ai++] = basis[2].x;
					outArray[ai++] = basis[2].y;
					outArray[ai++] = basis[2].z;

					outArray[ai++] = p.x;
					outArray[ai++] = p.y;
					outArray[ai++] = p.z;

					outArray[ai++] = sphere.center.x;
					outArray[ai++] = sphere.center.y;
					outArray[ai++] = sphere.center.z;
					outArray[ai++] = sphere.radius;
				}
				
				// add plane
				outArray[ai++] = x;
				outArray[ai++] = y;
				outArray[ai++] = theta;
				outArray[ai++] = phi;

				outArray[ai++] = basis[1].x;
				outArray[ai++] = basis[1].y;
				outArray[ai++] = basis[1].z;

				outArray[ai++] = basis[2].x;
				outArray[ai++] = basis[2].y;
				outArray[ai++] = basis[2].z;

				outArray[ai++] = p.x;
				outArray[ai++] = p.y;
				outArray[ai++] = p.z;

				outArray[ai++] = 0;
				outArray[ai++] = 0;
				outArray[ai++] = 0;
				outArray[ai++] = -1;

/*
				// Select a random direction
				var x = cos(phi) * theta;
				var y = sin(phi) * theta;
				var z = sqrt(1.0 - theta * theta);
				
				// Transform ray direction on local plane to global coordinate
				var rx = x * basis[0].x + y * basis[1].x + z * basis[2].x;
				var ry = x * basis[0].y + y * basis[1].y + z * basis[2].y;
				var rz = x * basis[0].z + y * basis[1].z + z * basis[2].z;
*/				
			}
		}
		
		return ai;
	}

	function parCalcAO(index) {
		var topIndex = index[0] - 0;
		var el = this.get(topIndex);

		var sx = el[0];
		var sy = el[1];
		var xy = sx*1000 + sy;

		var theta = el[2];
		var phi = el[3];

		var x = Math.cos(phi) * theta;
		var y = Math.sin(phi) * theta;
		var z = Math.sqrt(1.0 - theta * theta);
		
		var bY_x = el[4];
		var bY_y = el[5];
		var bY_z = el[6];

		var bZ_x = el[7];
		var bZ_y = el[8];
		var bZ_z = el[9];
		
		var ox = el[10];
		var oy = el[11];
		var oz = el[12];

		var cx = el[13];
		var cy = el[14];
		var cz = el[15];
		var R  = el[16];
		var occ = 0;
		
		// bX <- bY X bZ
		var bX_x = bY_y * bZ_z - bY_z * bZ_y;
		var bX_y = bY_z * bZ_x - bY_x * bZ_z;
		var bX_z = bY_x * bZ_y - bY_y * bZ_x;

		var rdx = x * bX_x + y * bY_x + z * bZ_x;
		var rdy = x * bX_y + y * bY_y + z * bZ_y;
		var rdz = x * bX_z + y * bY_z + z * bZ_z;

		if (R < 0) {
			// Plane
			var ppx = 0.0;
			var ppy = -0.5;
			var ppz = 0.0;

			var pnx = 0.0;
			var pny = 1.0;
			var pnz = 0.0;
			
			var pd = -(ppx * pnx + ppy * pny + ppz * pnz);
			var pv = rdx * pnx + rdy * pny + rdz * pnz;
			var pva = (pv < 0) ? -pv : pv;

			if (pva >= 1.0e-17) {
				var pt = -((ox * pnx + oy * pny + oz * pnz) + pd) / pv;
				if (pt > 0.0) {
					occ = 1;
				}
			}
		} else {
			// Sphere
			
		
			/* sphere isect */
			var rsx = ox - cx;
			var rsy = oy - cy;
			var rsz = oz - cz;

			var B = rsx*rdx + rsy*rdy + rsz*rdz;
			var C = rsx*rsx + rsy*rsy + rsz*rsz - R * R;
			var D = B * B - C;

			if (D > 0.0) {
				var s_t = -B - Math.sqrt(D);	
				if (s_t > 0.0) {
					occ = 1;
				}
			}
		}
		return [sx, sy, occ];
	}

	function parRender(img, w, h, nsubsamples){
		
		
		
		// Work areas (to suppress mass allocation)
		var rsTempVec = new Vec();
		var aoRayTemp = new Ray();
		var aoOrgTemp = new Vec();
		var aoIsectTemp = new Isect();
		var basisTemp = [new Vec(), new Vec(), new Vec()];

		var col = new Vec();
	    var x, y;
	    var u, v;

		var ray = new Ray();
		var isect = new Isect();
	    var fimg = allocateFloats(w * h * 3);
		var aoQueue = new Float32Array(w*h * nsubsamples*nsubsamples * 4 * NAO_SAMPLES*NAO_SAMPLES * 17);
		var aqIndex = 0;

		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				for (v = 0; v < nsubsamples; ++v) {
					for (u = 0; u < nsubsamples; ++u) {
						// Make normalized coordinate
						var px = (x + (u / nsubsamples) - (w / 2.0)) / (w / 2.0);
						var py = -(y + (v / nsubsamples) - (h / 2.0)) / (h / 2.0);
						
						// Eye position
						ray.org.x = 0.0;
						ray.org.y = 0.0;
						ray.org.z = 0.0;
						
						// Ray direction from the eye
						ray.dir.x = px;
						ray.dir.y = py;
						ray.dir.z = -1.0;
						ray.dir.normalize();
						
						isect.t   = 1.0e+17;
						isect.hit = false;
						
						// Cast the primary ray
						raySphereIntersect(isect, ray, gScene.spheres[0], rsTempVec);
						raySphereIntersect(isect, ray, gScene.spheres[1], rsTempVec);
						raySphereIntersect(isect, ray, gScene.spheres[2], rsTempVec);
						rayPlaneIntersect (isect, ray, gScene.plane);

						if (isect.hit) {
							// Do secondary ray tracing
//							ambientOcclusion(gScene, col, isect, rsTempVec, aoRayTemp, aoOrgTemp, aoIsectTemp, basisTemp);
							aqIndex = addAOEntry(aoQueue, aqIndex, gScene, x, y, isect, aoOrgTemp, basisTemp);
//							aoQueue.push(x,y, isect.p.x, isect.p.y, isect.p.z, isect.n.x, isect.n.y, isect.n.z);
							/*
							fimg[3 * (y * w + x) + 0] += col.x;
							fimg[3 * (y * w + x) + 1] += col.y;
							fimg[3 * (y * w + x) + 2] += col.z;
							*/
						}
					}
				}
				/*
				var bufpos = 3 * (y * w + x);
				fimg[bufpos + 0] /= nsubsamples * nsubsamples;
				fimg[bufpos + 1] /= nsubsamples * nsubsamples;
				fimg[bufpos + 2] /= nsubsamples * nsubsamples;
				
				img[bufpos + 0] = fimg[bufpos + 0] * 255;
				img[bufpos + 1] = fimg[bufpos + 1] * 255;
				img[bufpos + 2] = fimg[bufpos + 2] * 255;
				*/
			}
		}
	
		var paAO = new ParallelArray(aoQueue);
//		paAO.shape = [paAO.length/17, 17];
//		paAO.strides = [17, 1];
		var paAOResult  = paAO.partition(17).combine(parCalcAO);
		var aoLenPerSample = (NAO_SAMPLES*NAO_SAMPLES) * 4; 
		var paOccCounts = paAOResult.partition(aoLenPerSample).combine(function(index) {
			var a = this.get(index[0] - 0);
			var nDirs = a.length >> 2;
			var occ = 0;
			for (var di = 0;di < nDirs;di++) {
				var oi = di << 2;
				if (a[oi][2] | a[oi+1][2] > 0 | a[oi+2][2] > 0 | a[oi+3][2]) {
					occ++;
				}
			}
			
			var iocc = (nDirs - occ) / nDirs;
			return [a[0][0], a[0][1], iocc];
		});
		
		
		var q;
		/*
		var cnt = 0;
		for (q = 0;q < 190;++q) {
			var occ = paOccCounts.get(30000+q);
			console.log(occ);
		}
		*/
		
		var oIndex = 0;
		var intensity;
//		console.log(paOccCounts)
//		return;
		
		var farr = paOccCounts.flatten();
		for (q = 0;q < 16384;++q) {
			for (v = 0; v < nsubsamples; ++v) {
				for (u = 0; u < nsubsamples; ++u) {
					x = farr.get(oIndex++) | 0;
					y = farr.get(oIndex++) | 0;
					var intensity = farr.get(oIndex++);
					fimg[3 * (y * w + x) + 0] += intensity;
					fimg[3 * (y * w + x) + 1] += intensity;
					fimg[3 * (y * w + x) + 2] += intensity;
				}
			}
//			console.log(x, y);
			/*
*/
		}
		
		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				var bufpos = 3 * (y * w + x);
				fimg[bufpos + 0] /= nsubsamples * nsubsamples;
				fimg[bufpos + 1] /= nsubsamples * nsubsamples;
				fimg[bufpos + 2] /= nsubsamples * nsubsamples;
				
				img[bufpos + 0] = fimg[bufpos + 0] * 255;
				img[bufpos + 1] = fimg[bufpos + 1] * 255;
				img[bufpos + 2] = fimg[bufpos + 2] * 255;
			}
		}
	}


	/*
	
	
	function parRender(img, w, h, nsubsamples)
	{
		var spheres = gScene.spheres;
		var N_SPHERES = spheres.length;
		var RENDER_ELEMENT_SIZE = 12;
		var x, y, i;
		var u, v, s;
		var isect = new Isect();
		var paIndex = [0, 0];
		var col = new Vec();
		var ray = new Ray();

	    var fimg = allocateFloats(w * h * 3);
		var elem_count = w * h * nsubsamples * nsubsamples * N_SPHERES;
		var paSource = new Float32Array(elem_count * RENDER_ELEMENT_SIZE);

		// Eye position
		ray.org.x = 0.0;
		ray.org.y = 0.0;
		ray.org.z = 0.0;

		var paPos = 0;
		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				for (v = 0; v < nsubsamples; ++v) {
	                for (u = 0; u < nsubsamples; ++u) {
						for (s = 0;s < N_SPHERES;++s) {

							// Make normalized coordinate
							var px = (x + (u / nsubsamples) - (w / 2.0)) / (w / 2.0);
							var py = -(y + (v / nsubsamples) - (h / 2.0)) / (h / 2.0);

							// Ray direction from the eye
							ray.dir.x = px;
							ray.dir.y = py;
							ray.dir.z = -1.0;
							ray.dir.normalize();
						
							// Eye position
							paSource[paPos++] = 0;
							paSource[paPos++] = 0;
							paSource[paPos++] = 0;
							paSource[paPos++] = 0;
							
							// Ray direction from the eye
							paSource[paPos++] = ray.dir.x;
							paSource[paPos++] = ray.dir.y;
							paSource[paPos++] = ray.dir.z;
							paSource[paPos++] = 0;

							paSource[paPos++] = spheres[s].center.x;
							paSource[paPos++] = spheres[s].center.y;
							paSource[paPos++] = spheres[s].center.z;
							paSource[paPos++] = spheres[s].radius;
						}
					}
				}
				
			} // end line
		}
		var paSphereRays = new ParallelArray(paSource);
		paSphereRays.shape = [elem_count, RENDER_ELEMENT_SIZE];
		paSphereRays.strides = [RENDER_ELEMENT_SIZE, 1];
		var paSphereRaysResult  = paSphereRays.combine(parRaytraceSphere);
		var flatRes = paSphereRaysResult.flatten();
		return
///
		var RESULT_ELEMENT_SIZE = 7;
		paPos = 0;
		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				for (v = 0; v < nsubsamples; ++v) {
					for (u = 0; u < nsubsamples; ++u) {
						isect.t = 1.0e+17;
						isect.hit = false;
						col.x = col.y = col.z = 0;

						for (s = 0;s < N_SPHERES;++s) {
							var st = flatRes.get(paPos);
//							console.log(st)
							if (st > 0 && st < isect.t) {
								isect.t = st;
								isect.hit = true;
								
								isect.p.x = flatRes.get(paPos+1);
								isect.p.y = flatRes.get(paPos+2);
								isect.p.z = flatRes.get(paPos+3);

								isect.n.x = flatRes.get(paPos+4);
								isect.n.y = flatRes.get(paPos+5);
								isect.n.z = flatRes.get(paPos+6);
							}
							paPos += RESULT_ELEMENT_SIZE;
						}


						// Ray direction from the eye
						ray.dir.x = (x + (u / nsubsamples) - (w / 2.0)) / (w / 2.0);
						ray.dir.y = -(y + (v / nsubsamples) - (h / 2.0)) / (h / 2.0);
						ray.dir.z = -1.0;
						ray.dir.normalize();
						rayPlaneIntersect (isect, ray, gScene.plane);

						if (isect.hit) {
							col.x = 1;
							col.y = 1;
							col.z = 1;
						}

						fimg[3 * (y * w + x) + 0] += col.x;
						fimg[3 * (y * w + x) + 1] += col.y;
						fimg[3 * (y * w + x) + 2] += col.z;

					}
				}
				
				var bufpos = 3 * (y * w + x);
				fimg[bufpos + 0] /= nsubsamples * nsubsamples;
				fimg[bufpos + 1] /= nsubsamples * nsubsamples;
				fimg[bufpos + 2] /= nsubsamples * nsubsamples;
				
				img[bufpos + 0] = fimg[bufpos + 0] * 255;
				img[bufpos + 1] = fimg[bufpos + 1] * 255;
				img[bufpos + 2] = fimg[bufpos + 2] * 255;

			} // end line
		}
		

	}
	*/
	
	function parRaytraceSphere(i) {
		var topIndex = i[0] - 0;
		
		// ----------------
		// fetch parameters
		var el = this.get(topIndex);
		var ox = el[0];
		var oy = el[1];
		var oz = el[2];

		var rdx = el[4];
		var rdy = el[5];
		var rdz = el[6];

		var cx = el[8];
		var cy = el[9];
		var cz = el[10];
		var R  = el[11];
		// ----------------
		
		var rsx = ox - cx;
		var rsy = oy - cy;
		var rsz = oz - cz;
		
		var B = rsx*rdx + rsy*rdy + rsz*rdz;
		var C = rsx*rsx + rsy*rsy + rsz*rsz - R * R;
		var D = B * B - C;

		var ret = [-1, 0,0,0, 0,0,0]
		if (D > 0.0) {
			var t = -B - Math.sqrt(D);
			
			if (t > 0.0) {
				ret[0] = t;

				ret[1] = ox + rdx * t;
				ret[2] = oy + rdy * t;
				ret[3] = oz + rdz * t;

				ret[4] = ret[1] - cx;
				ret[5] = ret[2] - cy;
				ret[6] = ret[3] - cz;
				
				var nlen = Math.sqrt(ret[4]*ret[4] + ret[5]*ret[5] + ret[6]*ret[6]);
				ret[4] /= nlen;
				ret[5] /= nlen;
				ret[6] /= nlen;
			}
		}

		return ret;
	}
	
	
})();