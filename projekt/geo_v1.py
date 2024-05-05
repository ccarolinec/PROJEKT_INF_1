from math import sin, cos, sqrt, tan, atan, atan2, degrees, radians, pi
import sys
import numpy as np

o = object()

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            
            
    def plh2xyz(self, phi, lam, h):
        phi = radians(phi)
        lam = radians(lam)
        N = self.a/(1 - self.ecc2 * (sin(phi))**2)**(1/2)
        X = (N + h) * cos(phi) * cos(lam)
        Y = (N + h) * cos(phi) * sin(lam)
        Z = (N * (1 - self.ecc2) + h) * sin(phi)
        return(X, Y, Z)
    
    def xyz2neu (self, X, Y, Z):
        phi, lam , h = xyz2plh(x, y, z)
        R = np.array([[-np.sin(phi)*np.cos(lam), -np.sin(lam), np.cos(phi)*np.cos(lam)], 
                      [-np.sin(phi)*np.sin(lam), np.cos(lam), np.cos(phi)*np.sin(lam)], 
                      [np.cos(phi), 0, np.sin(phi)]])
        XYZ = X, Y, Z
        X0 = 1
        Y0 = 2
        Z0 = 3 #robocze
        W0 = X0, Y0, Z0
        wektor_odl = XYZ - W0
        dx = R.T @ wektor_odl
        n = dx[0]
        e = dx[1]
        u = dx[-1]
        return(n, e, u)
    
    def pl22000(self, phi, lam):
        
        if lam >= np.deg2rad(13.5) and lam < np.deg2rad(16.5): #stopnie
            lam0 = np.deg2rad(15)
            nr = 5
        elif lam >=np.deg2rad(16.5) and lam < np.deg2rad(19.5):
            lam0 = np.deg2rad(18)
            nr = 6
        elif lam >=np.deg2rad(19.5) and lam < np.deg2rad(22.5):
            lam0 = np.deg2rad(21)
            nr = 7
        elif lam >=np.deg2rad(22.5) and lam < np.deg2rad(25.5):
            lam0 = np.deg2rad(24)
            nr = 8
        
         
        b2=self.a**2*(1-self.ecc2)
        e_2=(self.a**2-b2)/b2
        dl=lam-lam0
        t=tan(phi)
        n2=e_2*cos(phi)**2
        N=(self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2))
        A0=1-self.ecc2/4-3*self.ecc2**2/64-5*self.ecc2**3/256
        A2=(3/8)*(self.ecc2+self.ecc2**2/4+15*self.ecc2**3/128)
        A4=15/256*(self.ecc2**2+3*self.ecc2**3/4)
        A6=35*self.ecc2**3/3072
        sigma=self.a*(A0*phi-A2*sin(2*phi)+A4*sin(4*phi)-A6*sin(6*phi))
        x=sigma+(dl**2/2)*N*sin(phi)*cos(phi)*(1+(dl**2/12)*(cos(phi))**2*(5-t**2+9*n2+4*n2**2)+(dl**4/360)*(cos(phi))**4*(61-58*t**2+t**4+270*n2-330*n2*t**2))
        y=dl*N*cos(phi)*(1+(dl**2/6)*(cos(phi))**2*(1-t**2+n2)+(dl**4/120)*(cos(phi))**4*(5-18*t**2+t**4+14*n2-58*n2*t**2))
    
        x2000=x*0.999923
        y2000=y*0.999923+((nr*1000000)+500000)    
        
        return(x2000, y2000)
    
    def pl21992(self, phi, lam):
        lam0 = np.deg2rad(19)
            
        b2=self.a**2*(1-self.ecc2)
        e_2=(self.a**2-b2)/b2
        dl=lam-lam0
        t=tan(phi)
        n2=e_2*cos(phi)**2
        N=(self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2))
        A0=1-self.ecc2/4-3*self.ecc2**2/64-5*self.ecc2**3/256
        A2=(3/8)*(self.ecc2+self.ecc2**2/4+15*self.ecc2**3/128)
        A4=15/256*(self.ecc2**2+3*self.ecc2**3/4)
        A6=35*self.ecc2**3/3072
        sigma=self.a*(A0*phi-A2*sin(2*phi)+A4*sin(4*phi)-A6*sin(6*phi))
        x=sigma+(dl**2/2)*N*sin(phi)*cos(phi)*(1+(dl**2/12)*(cos(phi))**2*(5-t**2+9*n2+4*n2**2)+(dl**4/360)*(cos(phi))**4*(61-58*t**2+t**4+270*n2-330*n2*t**2))
        y=dl*N*cos(phi)*(1+(dl**2/6)*(cos(phi))**2*(1-t**2+n2)+(dl**4/120)*(cos(phi))**4*(5-18*t**2+t**4+14*n2-58*n2*t**2))
        
        m1992=0.9993
        x1992=x*m1992-5300000
        y1992=y*m1992+500000
        
        return(x1992, y1992)


if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    # X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    # phi, lam, h = geo.xyz2plh(X, Y, Z)
    # print(phi, lam, h)
    # phi, lam, h = geo.xyz2plh2(X, Y, Z)
    # print(phi, lam, h)
    print(sys.argv)
    input_file_path = sys.argv[-1]
    if '--xyz2plh' in sys.argv and '--phl2xyz' in sys.argv:
        print('możesz podać tylko jedną falgę')
    elif '--xyz2plh' in sys.argv:
            #1
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[4:]
            #print(coords_lines)
            
            coords_plh = []
            
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                phi, lam, h = geo.xyz2plh(x, y, z)
                coords_plh.append([phi, lam, h])
            
            
        with open('wsp_inp/result_xyz2plh.txt', 'w') as f:
            f.write('phi[deg], lam[deg], h[m]\n')
            
            for coords_list in coords_plh:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')
        
    elif '--plh2xyz' in sys.argv:
        
    #2
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[1:]
            #print(coords_lines)
            
            coords_xyz = []
            
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str, h_str = coord_line.split(',')
                phi, lam, h = (float(phi_str,), float(lam_str), float(h_str))
                x, y, z = geo.plh2xyz(phi, lam, h)
                coords_plh.append([x, y, z])
            
            
        with open('wsp_inp/result_plh2xyz.txt', 'w') as f:
            f.write('x[m], y[m], z[m]\n')
            
            for coords_list in coords_plh:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')
        

        
        
        
    
