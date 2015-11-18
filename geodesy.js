/* Contains code from OpenLayers project.
*  Copyright (c) 2006-2013 by OpenLayers Contributors (see authors.txt for
 * full list of contributors). Published under the 2-clause BSD license.
 * See license.txt in the OpenLayers distribution or repository for the
 * full text of the license. */

function Geodesy() {
    var self = this;
    self.rad =  function(x) {return x*Math.PI/180;};
    self.vincentyConstants = {
        a: 6378137,
        b: 6356752.3142,
        f: 1/298.257223563
    };
    self.distVincenty = function(p1, p2) {
        var ct = self.vincentyConstants;
        var a = ct.a, b = ct.b, f = ct.f;
        var L = self.rad(p2.lon - p1.lon);
        var U1 = Math.atan((1-f) * Math.tan(self.rad(p1.lat)));
        var U2 = Math.atan((1-f) * Math.tan(self.rad(p2.lat)));
        var sinU1 = Math.sin(U1), cosU1 = Math.cos(U1);
        var sinU2 = Math.sin(U2), cosU2 = Math.cos(U2);
        var lambda = L, lambdaP = 2*Math.PI;
        var iterLimit = 20;
        while (Math.abs(lambda-lambdaP) > 1e-12 && --iterLimit>0) {
            var sinLambda = Math.sin(lambda), cosLambda = Math.cos(lambda);
            var sinSigma = Math.sqrt((cosU2*sinLambda) * (cosU2*sinLambda) +
            (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda));
            if (sinSigma==0) {
                return 0;  // co-incident points
            }
            var cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
            var sigma = Math.atan2(sinSigma, cosSigma);
            var alpha = Math.asin(cosU1 * cosU2 * sinLambda / sinSigma);
            var cosSqAlpha = Math.cos(alpha) * Math.cos(alpha);
            var cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha;
            var C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
            lambdaP = lambda;
            lambda = L + (1-C) * f * Math.sin(alpha) *
            (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
        }
        if (iterLimit==0) {
            return NaN;  // formula failed to converge
        }
        var uSq = cosSqAlpha * (a*a - b*b) / (b*b);
        var A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
        var B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
        var deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
            B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
        var s = b*A*(sigma-deltaSigma);
        var d = s.toFixed(3)/1000; // round to 1mm precision
        return d;
    }
    /*
    Calculates geodesic length of a LineString.
    @param vertices Array of LatLng objects (with properties "lat" and "lng")
    @return float Length in meters
    */
    self.lineLength = function(vertices) {
        var length = 0.0;
        if(vertices.length > 1) {
            var p1, p2;
            for(var i=1, len=vertices.length; i<len; i++) {
                p1 = vertices[i-1];
                p2 = vertices[i];
                // this returns km and requires lon/lat properties
                length += self.distVincenty(
                    {lon: p1.lng, lat: p1.lat}, {lon: p2.lng, lat: p2.lat}
                );
            }
        }
        // convert to m
        return length * 1000;
    };
    /*
    Calculates geodesic area of a LinearRing.
    @param vertices Array of LatLng objects (with properties "lat" and "lng")
    @return float Area in sqare meters.
    */
    self.polygonArea = function(vertices) {
        var area = 0.0;
        var len = vertices.length;
        if(len > 2) {
            var p1, p2;
            for(var i=0; i<len-1; i++) {
                p1 = vertices[i];
                p2 = vertices[i+1];
                area += self.rad(p2.lng - p1.lng) *
                        (2 + Math.sin(self.rad(p1.lat)) +
                        Math.sin(self.rad(p2.lat)));
            }
            area = area * 6378137.0 * 6378137.0 / 2.0;
        }
        //the calculated area might be negative, so return absolute value
        return Math.abs(area);
    }
}
