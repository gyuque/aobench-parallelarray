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
	
	var WIDTH  = 192;
	var HEIGHT = 192;
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
	function packVector(v) {
		vx = ((v.x * 499) + 500) | 0;
		vy = ((v.y * 499) + 500) | 0;
		vz = ((v.z * 499) + 500) | 0;
		
		return vx + (vy * 1000) + (vz * 1000000)
	}
	
	function addAOEntry(outArray, ai, scn, x, y, isect, aoOrgTemp) {
		var i, j;
		var ntheta = NAO_SAMPLES;
		var nphi   = NAO_SAMPLES;
		var eps = 0.0001;

		var p = aoOrgTemp;
		p.x = isect.p.x + eps * isect.n.x;
		p.y = isect.p.y + eps * isect.n.y;
		p.z = isect.p.z + eps * isect.n.z;
/*		
		var basis = basisTemp;
		orthoBasis(basis, isect.n);
*/		
		var pi2 = M_PI * 2.0;

		// var sphereList = scn.spheres;
		// Do monte carlo sampling for secondary rays
		for (j = 0; j < ntheta; ++j) {
			for (i = 0; i < nphi; ++i) {
				var theta = sqrt(drand48());
				var phi   = pi2 * drand48();

				outArray[ai++] = x | (y << 12);
				outArray[ai++] = theta;
				outArray[ai++] = phi;

				outArray[ai++] = isect.n.x;
				outArray[ai++] = isect.n.y;
				outArray[ai++] = isect.n.z;

				outArray[ai++] = p.x;
				outArray[ai++] = p.y;
				outArray[ai++] = p.z;

				/*
				for (var si = 0;si < 3;++si) {
					var sphere = sphereList[si];
					
					outArray[ai++] = x;
					outArray[ai++] = y;
					outArray[ai++] = theta;
					outArray[ai++] = phi;

					outArray[ai++] = isect.n.x;
					outArray[ai++] = isect.n.y;
					outArray[ai++] = isect.n.z;

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

				outArray[ai++] = isect.n.x;
				outArray[ai++] = isect.n.y;
				outArray[ai++] = isect.n.z;

				outArray[ai++] = p.x;
				outArray[ai++] = p.y;
				outArray[ai++] = p.z;

				outArray[ai++] = 0;
				outArray[ai++] = 0;
				outArray[ai++] = 0;
				outArray[ai++] = -1;
				*/
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

		var sxsy = el[0];
		if (sxsy == 0) {return[0,0];}

		var theta = el[1];
		var phi = el[2];

		var x = Math.cos(phi) * theta;
		var y = Math.sin(phi) * theta;
		var z = Math.sqrt(1.0 - theta * theta);
		
		var bZ_x = el[3];
		var bZ_y = el[4];
		var bZ_z = el[5];

		var bY_x = 0;
		var bY_y = 0;
		var bY_z = 0;
		if ((bZ_x < 0.6) && (bZ_x > -0.6)) {
			bY_x = 1.0;
		} else if ((bZ_y < 0.6) && (bZ_y > -0.6)) {
			bY_y = 1.0;
		} else if ((bZ_z < 0.6) && (bZ_z > -0.6)) {
			bY_z = 1.0;
		} else {
			bY_x = 1.0;
		}

		var bX_x = bY_y * bZ_z - bY_z * bZ_y;
		var bX_y = bY_z * bZ_x - bY_x * bZ_z;
		var bX_z = bY_x * bZ_y - bY_y * bZ_x;
		var bXlen = Math.sqrt(bX_x*bX_x + bX_y*bX_y + bX_z*bX_z);
		bX_x /= bXlen;
		bX_y /= bXlen;
		bX_z /= bXlen;

		var x2 = bZ_y * bX_z - bZ_z * bX_y;
		var y2 = bZ_z * bX_x - bZ_x * bX_z;
		var z2 = bZ_x * bX_y - bZ_y * bX_x;
		var bYlen = Math.sqrt(x2*x2 + y2*y2 + z2*z2);
		
		x2 /= bYlen;
		y2 /= bYlen;
		z2 /= bYlen;
		
		
		var ox = el[6];
		var oy = el[7];
		var oz = el[8];

		var rdx = x * bX_x + y * x2 + z * bZ_x;
		var rdy = x * bX_y + y * y2 + z * bZ_y;
		var rdz = x * bX_z + y * z2 + z * bZ_z;
		var occ = 0;

		var R  = 0.5;
		for (var objIndex = 0;objIndex < 4 && !occ;objIndex++) {
			if (objIndex > 2) {
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
/*			
				var cx = -3.0 + objIndex * 1.5;
				var cy = 0;
				var cz = (objIndex == 0) ? -3.5 :
				         (objIndex == 1) ? -3.0 : 
				         2.2;
				var R  = 0.5;
*/
				var cx, cy, cz;
				if (objIndex == 0) {
					cx = -2.0;
					cy = 0;
					cz = -3.5;
				} else if (objIndex == 1) {
					cx = -0.5;
					cy = 0;
					cz = -3;
				} else {
					cx = 1.0;
					cy = 0;
					cz = -2.2;
				}

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
			
		}
		
		return [sxsy, occ];
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
		var aoQueue = new Float32Array(w*h * nsubsamples*nsubsamples * NAO_SAMPLES*NAO_SAMPLES * 9);
		console.log(aoQueue.length)
		var aqIndex = 0;
		var aoCount = 0;

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
							// Prepare secondary ray tracing
							aqIndex = addAOEntry(aoQueue, aqIndex, gScene, x, y, isect, aoOrgTemp);
							++aoCount;
						}
					}
				}

			}
		}

		var aoLenPerSample = NAO_SAMPLES*NAO_SAMPLES; 
		var paAO = new ParallelArray(aoQueue);
		var paAOResult  = paAO.partition(9).combine(parCalcAO);

		var paOccCounts = paAOResult.partition(aoLenPerSample).combine(function(index) {
			var a = this.get(index[0] - 0);
			var nDirs = a.length;
			var occ = 0;
			for (var di = 0;di < nDirs;di++) {
				occ += a[di][1];
			}
			
			var iocc = (nDirs - occ) / nDirs;
			return [a[0][0], iocc];
		});

		var q;
		var oIndex = 0;
		var intensity;
		var farr = paOccCounts.flatten();
		
		for (q = 0;q < aoCount;++q) {
			var xy = farr.get(oIndex++) | 0;
			x = xy & 0xfff;
			y = xy >> 12;
			var intensity = farr.get(oIndex++);
			fimg[3 * (y * w + x) + 0] += intensity;
			fimg[3 * (y * w + x) + 1] += intensity;
			fimg[3 * (y * w + x) + 2] += intensity;
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

})();