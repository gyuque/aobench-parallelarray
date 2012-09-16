(function() {
	//'use strict';
	// ++++ Original version ++++
	
	// Function shortcuts
	var sqrt = Math.sqrt;
	var fabs = Math.abs;
	var drand48 = Math.random;
	var M_PI = Math.PI;
	var cos = Math.cos;
	var sin = Math.sin;
	
	var WIDTH  = 256;
	var HEIGHT = 256;
	var NSUBSAMPLES = 2;
	var NAO_SAMPLES = 8;
	
	var gScene = null;
	
	
	// Geometries ---------------------------------
	/*
	function Vec() {
		this.x = 0;
		this.y = 0;
		this.z = 0;
	}
	*/
	function createVec() {
		return [0, 0, 0,  0];
	}
	
	function vdot(a, b) {
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}
	
	function vnormalize(v) {
		var length = sqrt(vdot(v, v));

		if (fabs(length) > 1.0e-17) {
			v[0] /= length;
			v[1] /= length;
			v[2] /= length;
		}
	}
	
	function vcross(vOut, v0, v1) {
		vOut[0] = v0[1] * v1[2] - v0[2] * v1[1];
		vOut[1] = v0[2] * v1[0] - v0[0] * v1[2];
		vOut[2] = v0[0] * v1[1] - v0[1] * v1[0];
	}

/*	
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
	*/
	
	function createSphere(x, y, z, r) {
		return [x, y, z, r];
	}
	
	// Rendering informations -----------------------
	function creataRay() {
		return [
			createVec(), // origin position
			createVec()  // normal
		]
	}
	
	function createIsect() {
		return [
			createVec(), // position
			createVec(), // normal
			createVec()  // distance, flag
		];
	}
	/*
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
	*/
	
	// Hit tests
	function raySphereIntersect(isectPos, isectN, isectF, rayPos, rayDir, sphere, tempVec) {
		var sphereR = sphere[3];
		var rs = tempVec;
		
		rs[0] = rayPos[0] - sphere[0];
		rs[1] = rayPos[1] - sphere[1];
		rs[2] = rayPos[2] - sphere[2];
		var B = vdot(rs, rayDir);
		var C = vdot(rs, rs) - sphereR * sphereR;
		var D = B * B - C;

		if (D > 0.0) {
			var t = -B - sqrt(D);
			if ((t > 0.0) && (t < isectF[0])) {
				isectF[0] = t;
				isectF[1] = 1;

				isectPos[0] = rayPos[0] + rayDir[0] * t;
				isectPos[1] = rayPos[1] + rayDir[1] * t;
				isectPos[2] = rayPos[2] + rayDir[2] * t;

				isectN[0] = isectPos[0] - sphere[0];
				isectN[1] = isectPos[1] - sphere[1];
				isectN[2] = isectPos[2] - sphere[2];

				vnormalize(isectN);
			}
		}
	}
	
	function rayPlaneIntersect(isectPos, isectN, isectF, rayPos, rayDir, planePos, planeN)
	{
		var d = -vdot(planePos, planeN);
		var v = vdot(rayDir, planeN);
//console.log()

		if (fabs(v) < 1.0e-17) return;

		var t = -(vdot(rayPos, planeN) + d) / v;
		if ((t > 0.0) && (t < isectF[0])) {
			isectF[0] = t;
			isectF[1] = 1;

			isectPos[0] = rayPos[0] + rayDir[0] * t;
			isectPos[1] = rayPos[1] + rayDir[1] * t;
			isectPos[2] = rayPos[2] + rayDir[2] * t;

			isectN[0] = planeN[0];
			isectN[1] = planeN[1];
			isectN[2] = planeN[2];
		}
	}

	// Renderer functions

	function orthoBasis(basisX, basisY, basisZ, n)
	{
		basisZ[0] = n[0];
		basisZ[1] = n[1];
		basisZ[2] = n[2];
		
		basisY[0] = 0.0; basisY[1] = 0.0; basisY[2] = 0.0;

		if ((n[0] < 0.6) && (n[0] > -0.6)) {
			basisY[0] = 1.0;
	    } else if ((n[1] < 0.6) && (n[1] > -0.6)) {
			basisY[1] = 1.0;
	    } else if ((n[2] < 0.6) && (n[2] > -0.6)) {
			basisY[2] = 1.0;
	    } else {
			basisY[0] = 1.0;
	    }

		vcross(basisX, basisY, basisZ);
		vnormalize(basisX);

		vcross(basisY, basisZ, basisX);
		vnormalize(basisY);
	}

	function ambientOcclusion(
			sceneArray, sceneStart,
			col, isectPos, isectN, isectF,
			rsTempVec, tmpRayPos, tmpRayDir,
			aoOrgTemp, tmpIsectPos, tmpIsectDir, tmpIsectF,
			basisArray, basisStart) {
		var i, j;
		var ntheta = NAO_SAMPLES;
		var nphi   = NAO_SAMPLES;
		var eps = 0.0001;

		var p = aoOrgTemp;
		p[0] = isectPos[0] + eps * isectN[0];
		p[1] = isectPos[1] + eps * isectN[1];
		p[2] = isectPos[2] + eps * isectN[2];
		
		var basisX = basisArray[basisStart    ];
		var basisY = basisArray[basisStart + 1];
		var basisZ = basisArray[basisStart + 2];
		orthoBasis(basisX, basisY, basisZ, isectN);
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
				var rx = x * basisX[0] + y * basisY[0] + z * basisZ[0];
				var ry = x * basisX[1] + y * basisY[1] + z * basisZ[1];
				var rz = x * basisX[2] + y * basisY[2] + z * basisZ[2];
	
				tmpRayPos[0] = p[0];
				tmpRayPos[1] = p[1];
				tmpRayPos[2] = p[2];
				tmpRayDir[0] = rx;
				tmpRayDir[1] = ry;
				tmpRayDir[2] = rz;
				
				tmpIsectF[0] = 1.0e+17;
				tmpIsectF[1] = 0;

				// Cast a ray
				raySphereIntersect(tmpIsectPos, tmpIsectDir, tmpIsectF, tmpRayPos, tmpRayDir, sceneArray[sceneStart  ], rsTempVec);
				raySphereIntersect(tmpIsectPos, tmpIsectDir, tmpIsectF, tmpRayPos, tmpRayDir, sceneArray[sceneStart+1], rsTempVec);
				raySphereIntersect(tmpIsectPos, tmpIsectDir, tmpIsectF, tmpRayPos, tmpRayDir, sceneArray[sceneStart+2], rsTempVec);
				rayPlaneIntersect(tmpIsectPos, tmpIsectDir, tmpIsectF, tmpRayPos, tmpRayDir, sceneArray[sceneStart+3], sceneArray[sceneStart+4]);
				if (tmpIsectF[1]) occlusion += 1.0;
			}
		}
		
		occlusion = (ntheta * nphi - occlusion) / (ntheta * nphi);

		// Final irradiance
		col[0] = occlusion;
		col[1] = occlusion;
		col[2] = occlusion;
	}	
	
	function initScene()
	{
		var scene = [
			createSphere(-2.0, 0.0, -3.5, 0.5),
			createSphere(-0.5, 0.0, -3.0, 0.5),
			createSphere(1.0, 0.0, -2.2, 0.5),
			[0.0, -0.5, 0.0, 0],
			[0.0, 1.0, 0.0, 0]
		];
	
		return scene;
	}
	
	function render(img, w, h, nsubsamples)
	{
		// Work areas (to suppress mass allocation)
		var rsTempVec = createVec();
		var aoRayTemp = creataRay();
		var aoOrgTemp = createVec();
		var aoIsectTemp = createIsect();
		var basisTemp = [createVec(), createVec(), createVec()];

		var col = createVec();
	    var x, y;
	    var u, v;

		var ray = creataRay();
		var isect = createIsect();
		var fimg = allocateFloats(w * h * 3);
		var scn = gScene;

		for (y = 0; y < h; ++y) {
			for (x = 0; x < w; ++x) {
				for (v = 0; v < nsubsamples; ++v) {
	                for (u = 0; u < nsubsamples; ++u) {
						// Make normalized coordinate
						var px = (x + (u / nsubsamples) - (w / 2.0)) / (w / 2.0);
						var py = -(y + (v / nsubsamples) - (h / 2.0)) / (h / 2.0);
						
						// Eye position
						ray[0][0] = 0.0;
						ray[0][1] = 0.0;
						ray[0][2] = 0.0;
						
						// Ray direction from the eye
						ray[1][0] = px;
						ray[1][1] = py;
						ray[1][2] = -1.0;
						vnormalize(ray[1]);
						
						isect[2][0] = 1.0e+17;
						isect[2][1] = 0;
						
						// Cast the primary ray
						raySphereIntersect(isect[0], isect[1], isect[2], ray[0], ray[1], scn[0], rsTempVec);
						raySphereIntersect(isect[0], isect[1], isect[2], ray[0], ray[1], scn[1], rsTempVec);
						raySphereIntersect(isect[0], isect[1], isect[2], ray[0], ray[1], scn[2], rsTempVec);
						rayPlaneIntersect(isect[0], isect[1], isect[2], ray[0], ray[1], scn[3], scn[4]);

						if (isect[2][1]) {
							// Do secondary ray tracing
							ambientOcclusion(
								scn, 0,
								col,
								isect[0], isect[1], isect[2],
								rsTempVec, aoRayTemp[0], aoRayTemp[1],
								aoOrgTemp, aoIsectTemp[0], aoIsectTemp[1], aoIsectTemp[2],
								basisTemp, 0);
							fimg[3 * (y * w + x) + 0] += col[0];
							fimg[3 * (y * w + x) + 1] += col[1];
							fimg[3 * (y * w + x) + 2] += col[2];
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
	
	function parRender(img, w, h, nsubsamples)
	{
	    var x, y, i;
	    var u, v;

		var vec = createVec();
	    var fimg = allocateFloats(w * h * 3);
		function parRaytrace(elem) {
			var scene = []
			/*
			throw 0;
			var isect = elem.isect;
			var work = elem.workArea;
			var rsTempVec = work.rsVec;
			var scn = elem.scene;
			
			raySphereIntersect(isect, elem.ray, scn.spheres[0], rsTempVec);
			raySphereIntersect(isect, elem.ray, scn.spheres[1], rsTempVec);
			raySphereIntersect(isect, elem.ray, scn.spheres[2], rsTempVec);
			rayPlaneIntersect (isect, elem.ray, scn.plane);
			
			if (isect.hit) {
				ambientOcclusion(scn, elem.color, isect, rsTempVec, work.aoRay, work.aoOrg, work.aoIsect, work.aoBasis);
			} 
			*/
			return elem;
		}

		var paSource = makeRenderWorkArray(w * nsubsamples * nsubsamples, gScene);
		for (y = 0; y < h; ++y) {
			// Setup rays
			var paPos = 0;
			for (x = 0; x < w; ++x) {
				for (v = 0; v < nsubsamples; ++v) {
	                for (u = 0; u < nsubsamples; ++u) {
						// Make normalized coordinate
						var px = (x + (u / nsubsamples) - (w / 2.0)) / (w / 2.0);
						var py = -(y + (v / nsubsamples) - (h / 2.0)) / (h / 2.0);
						
						// Ray
						// Eye position
						paSource[paPos++] = 0.0;
						paSource[paPos++] = 0.0;
						paSource[paPos++] = 0.0;
						++paPos;
						
						// Ray direction from the eye
						vec[0] = px;
						vec[1] = py;
						vec[2] = -1.0;
						vnormalize(vec);
						paSource[paPos++] = vec[0];
						paSource[paPos++] = vec[1];
						paSource[paPos++] = vec[2];
						++paPos;

						// Isect
						paPos += 4;
						paPos += 4;
						
						paSource[paPos++] = 1.0e+17;
						paSource[paPos++] = 0;
						paPos += 2;

						// Out color
						paSource[paPos++] = 0;
						paSource[paPos++] = 0;
						paSource[paPos++] = 0;
						++paPos;
						
						paPos += 60;
					}
				}
				
			} // end line
//console.log(paPos, w * nsubsamples * nsubsamples * 84)
			var parArray = new ParallelArray(paSource);
			var paResults = parArray.map(parRaytrace);
			throw 456;
			
			var paPos = 0;
			for (x = 0; x < w; ++x) {
				var n2 = nsubsamples * nsubsamples;
				for (var i = 0;i < n2;++i) {
					var outCol = paResults[paPos++].color;
					fimg[3 * (y * w + x) + 0] += outCol.x;
					fimg[3 * (y * w + x) + 1] += outCol.y;
					fimg[3 * (y * w + x) + 2] += outCol.z;
				}

				var bufpos = 3 * (y * w + x);
				fimg[bufpos + 0] /= n2;
				fimg[bufpos + 1] /= n2;
				fimg[bufpos + 2] /= n2;
				
				img[bufpos + 0] = fimg[bufpos + 0] * 255;
				img[bufpos + 1] = fimg[bufpos + 1] * 255;
				img[bufpos + 2] = fimg[bufpos + 2] * 255;
			}
		}
	}
	
	function createPARenderElement(scene) {
		return [
			// Ray
			/* 0*/ 0,0,0,0, // origin position
			/* 1*/ 0,0,0,0, // normal
			
			// Isect
			/* 2*/ 0,0,0,0, // position
			/* 3*/ 0,0,0,0, // normal
			/* 4*/ 0,0,0,0, // distance, flag

			// Out color
			/* 5*/ 0,0,0,0,
			
			// Scene
			/* 6*/ scene[0][0], scene[0][1], scene[0][2], scene[0][3],
			/* 7*/ scene[1][0], scene[1][1], scene[1][2], scene[1][3],
			/* 8*/ scene[2][0], scene[2][1], scene[2][2], scene[2][3],
			/* 9*/ scene[3][0], scene[3][1], scene[3][2], scene[3][3],
			/*10*/ scene[4][0], scene[4][1], scene[4][2], scene[4][3],
			
			// Work area
			/*11*/ 0,0,0,0, // rsTempVec

			/*12*/ 0,0,0,0, // aoTempRay
			/*13*/ 0,0,0,0,
			
			/*14*/ 0,0,0,0, // aoTemoOrigin
			
			/*15*/ 0,0,0,0, // aoTempIsect
			/*16*/ 0,0,0,0,
			/*17*/ 0,0,0,0,
			
			/*18*/ 0,0,0,0, // aoTempBasis
			/*19*/ 0,0,0,0,
			/*20*/ 0,0,0,0
		];
	}
	
	function makeRenderWorkArray(length, scene) {
		// Setup elements
		var elements = [];
		for (var i = 0;i < length;++i) {
			Array.prototype.push.apply(elements, createPARenderElement(scene));
		}
		
		return elements;
		/*
		try {
			var pa = new ParallelArray(elements);
		} catch(e) {
			alert("Your browser does not support ParallelArray.");
			return null;
		}*/
		//return pa;
	}

})();