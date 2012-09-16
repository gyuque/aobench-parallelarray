(function() {
	'use strict';
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
	
	var _tempVec = new Vec();
	function raySphereIntersect(isect, ray, sphere) {
		var rs = _tempVec;
		
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

	var _aoRayTemp = new Ray();
	var _aoOrgTemp = new Vec();
	var _aoIsectTemp = new Isect();
	var _basisTemp = [new Vec(), new Vec(), new Vec()];
	function ambientOcclusion(col, isect) {
		var i, j;
		var ntheta = NAO_SAMPLES;
		var nphi   = NAO_SAMPLES;
		var eps = 0.0001;

		var p = _aoOrgTemp;
		p.x = isect.p.x + eps * isect.n.x;
		p.y = isect.p.y + eps * isect.n.y;
		p.z = isect.p.z + eps * isect.n.z;
		
		var basis = _basisTemp;
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
				
				var ray = _aoRayTemp;
				ray.org.x = p.x;
				ray.org.y = p.y;
				ray.org.z = p.z;
				ray.dir.x = rx;
				ray.dir.y = ry;
				ray.dir.z = rz;
				
				var occIsect = _aoIsectTemp;
				occIsect.t   = 1.0e+17;
				occIsect.hit = false;

				// Cast a ray
				raySphereIntersect(occIsect, ray, gScene.spheres[0]);
				raySphereIntersect(occIsect, ray, gScene.spheres[1]);
				raySphereIntersect(occIsect, ray, gScene.spheres[2]);
				rayPlaneIntersect (occIsect, ray, gScene.plane);
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
						raySphereIntersect(isect, ray, gScene.spheres[0]);
						raySphereIntersect(isect, ray, gScene.spheres[1]);
						raySphereIntersect(isect, ray, gScene.spheres[2]);
						rayPlaneIntersect (isect, ray, gScene.plane);

						if (isect.hit) {
							// Do secondary ray tracing
							ambientOcclusion(col, isect);

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
		function startAO() {
			var outCanvas = document.getElementById("out-canvas");
			var g = outCanvas.getContext("2d");
			outCanvas.width  = WIDTH;
			outCanvas.height = HEIGHT;
		
			var img = allocateBytes(WIDTH * HEIGHT * 3);
			gScene = initScene();
			render(img, WIDTH, HEIGHT, NSUBSAMPLES);
		
			emitToCanvas(g, img, WIDTH, HEIGHT);
		}
		
		var btn = document.getElementById("start-btn");
		btn.addEventListener('click', function() {
			btn.disabled = true;
			var startT = new Date();
			startAO();
			document.getElementById("result-out").innerHTML = (new Date() - startT) + " ms - " + 
			                                                  WIDTH+"x"+HEIGHT+" px, "+NSUBSAMPLES+"x"+NSUBSAMPLES+" subsamples";
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
						raySphereIntersect(isect, ray, gScene.spheres[0]);
						raySphereIntersect(isect, ray, gScene.spheres[1]);
						raySphereIntersect(isect, ray, gScene.spheres[2]);
						rayPlaneIntersect (isect, ray, gScene.plane);

						if (isect.hit) {
							// Do secondary ray tracing
							ambientOcclusion(col, isect);

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

})();