import svgwrite
import math

def rotate(x,y,angle,ox,oy):
    angle = math.radians(angle) #convert angle in deg to radians
    x_rot = ox + math.cos(angle) * (x - ox) - math.sin(angle) * (y - oy)
    y_rot = oy + math.sin(angle) * (x - ox) + math.cos(angle) * (y - oy)
    return x_rot,y_rot

def xy_angle(x1,y1,x2,y2):
    theta = (math.atan2(y2,x2) - math.atan2(y1,x1)) * (180.0 / math.pi)
    theta % 360
    return theta

def xy_distance(x1,y1,x2,y2):
    dist=math.sqrt(math.pow((x2-x1),2) + math.pow((y2-y1),2))
    return dist

class hcCanvas:
    x_size = 650
    y_size = 500
    axColor = "rgb(200,200,200)"
    axColorMin = "rgb(230,230,230)"
    axGonColor = "rgb(255,255,255)"
    axSW=2
    axSWGon=1
    axMax=20
    scalar = axMax*1/3+1; #to transform from pixels to units
    z_gap = (axMax+15) * scalar
    h_pentagon = (axMax * scalar * (math.sqrt(5+2*math.sqrt(5))/16))
    def __init__(self,svgfilename,hcType=''):
        self.filename = svgfilename
        self.hcType = hcType
        self.x_center = self.x_size/2-60
        self.y_center = self.y_size/2
        self.x_origin = self.x_center
        self.y_origin = self.y_size-self.y_center
        self.colors = options = {'uv':'rgb(127,15,126)','s':'rgb(11,36,251)','m':'rgb(15,127,18)','l':'rgb(252,13,27)','r':'rgb(0,0,0)','z':'rgb(253,189,64)'}
        self.canvas = svgwrite.Drawing(filename = self.filename + ".svg",size = (self.x_size, self.y_size), profile = "full", id = "eel")
        if hcType == 'H1': self.color_frame='rgb(255,101,128)'
        elif hcType == 'H2': self.color_frame='rgb(58,242,174)'
        elif hcType == 'H3': self.color_frame='rgb(255,216,0)'
        else: self.color_frame='rgb(0,0,0)'
        # self.canvas.add(self.canvas.rect(insert = (0,0),size = (self.x_size, self.y_size),id='frame',stroke="black",fill="rgb(240,240,240)"))
        self.canvas.add(self.canvas.rect(insert = (20,20),size = (self.x_size-20, self.y_size-20),stroke=self.color_frame,stroke_width=10,fill="white"))
#         self.addOrigin();
        self.addAxis();
        
        
    def addOrigin(self):
        self.origin=plotDatum(self,self.x_center,self.y_center,0,"none",self.axColor);
        self.origin.addDot('origin');
        self.z_origin=plotDatum(self,self.x_center+self.z_gap,self.y_center-(self.axMax * self.scalar),0,"none",self.axColor);
        self.z_origin.addDot('zOrigin');
    
    def drawAx(self,translation,tag=""): 
        if tag == 'uv': axLabel, axLabelAnchor, axLabelTrans = 'UV', 'middle', 'translate(-50,0)'
        elif tag == 's': axLabel, axLabelAnchor, axLabelTrans = 'S', 'start', 'translate(0,0)'
        elif tag == 'm': axLabel, axLabelAnchor, axLabelTrans = 'M', 'start', 'translate(50,0)'
        elif tag == 'l': axLabel, axLabelAnchor, axLabelTrans = 'L', 'end', 'translate(100,0)'
        elif tag == 'r': axLabel, axLabelAnchor, axLabelTrans = 'Rod', 'end', 'translate(150,0)'
        else: axLabel, axLabelAnchor, axLabelTrans = '', 'start', 'translate(200,0)'
        # draw axis    
        Ax = self.canvas.add(self.canvas.line(
            start=(self.x_origin,self.y_origin),
            end=(self.x_origin,self.y_origin-self.scalar*self.axMax),
            stroke=self.axColor,stroke_width=self.axSW,stroke_dasharray=4,
            transform='translate({0},0)'.format(translation)))
        # place label
        self.canvas.add(self.canvas.text(
            axLabel,
            insert=(rotate(self.x_origin+50, self.y_origin-self.scalar*self.axMax,0,self.x_origin,self.y_origin)),
            font_size=20,font_family='sans-serif',text_anchor=axLabelAnchor,
            fill=self.colors.get(tag), alignment_baseline="middle",
            transform=axLabelTrans
        ))
        return Ax
    
    def drawZAx(self):
        zAx = self.canvas.add(self.canvas.line(
            start=(self.x_origin,self.y_origin),
            end=(self.x_origin,self.y_origin-self.scalar*self.axMax),
            stroke=self.axColor,stroke_width=self.axSW,stroke_dasharray=4,
            transform='translate({0},0)'.format(250)))
        self.canvas.add(self.canvas.text(
            '?',
            insert=(rotate(self.x_origin+50, self.y_origin-self.scalar*self.axMax,0,self.x_origin,self.y_origin)),
            font_size=20,font_family='sans-serif',text_anchor='middle',
            fill=self.colors.get('z'),
            transform='translate(200,0)'
        ))
        
    def drawAxgon(self,scale):
        Axgon = self.canvas.add(
            svgwrite.shapes.Polygon(points=[
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*0/5,self.x_origin,self.y_origin),
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*1/5,self.x_origin,self.y_origin),
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*2/5,self.x_origin,self.y_origin),
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*3/5,self.x_origin,self.y_origin),
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*4/5,self.x_origin,self.y_origin),
            ],
                fill=self.axGonColor,stroke = self.axColor,stroke_width=self.axSWGon)
        )
        
    def drawAxgonMinor(self,scale):
        Axgon = self.canvas.add(
            svgwrite.shapes.Polygon(points=[
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*0/5,self.x_origin,self.y_origin),
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*1/5,self.x_origin,self.y_origin),
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*2/5,self.x_origin,self.y_origin),
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*3/5,self.x_origin,self.y_origin),
                rotate(self.x_origin,self.y_origin-self.scalar*scale,360*4/5,self.x_origin,self.y_origin),
            ],
                fill=self.axColorMin,stroke = self.axColor,stroke_width=self.axSWGon)
        )
    
    def drawZAxgon(self,scale):
        AxSquare = self.canvas.add(
            self.canvas.rect(
                insert = (self.x_center+self.z_gap+2.5-20,self.y_center+((self.axMax-scale) * self.scalar)-self.h_pentagon),
                size = (40,scale * self.scalar),
                fill=self.axGonColor,stroke = self.axColor,stroke_width=self.axSWGon)
        )
    
    def drawZAxgonMinor(self,scale):
        AxSquare = self.canvas.add(
            self.canvas.rect(
                insert = (self.x_center+self.z_gap+2.5-20,self.y_center+((self.axMax-scale) * self.scalar)-self.h_pentagon),
                size = (40,scale * self.scalar),
                fill=self.axColorMin,stroke = self.axColor,stroke_width=self.axSWGon)
        )
        
    
    
    def addAxis(self):
        # # self.drawAxgonMinor(25);
        # self.drawAxgonMinor(20);
        # self.drawAxgon(15);
        # self.drawAxgonMinor(10);
        # self.drawAxgon(5);
        # # self.drawZAxGonMinor(25);
        # self.drawZAxgonMinor(20);
        # self.drawZAxgon(15);
        # self.drawZAxgonMinor(10);
        # self.drawZAxgon(5);
        
#         [self.drawAx(360*facet/5) for facet in range(5)]
        self.drawZAx();
        self.axU=self.drawAx(50*0,'uv');
        self.axS=self.drawAx(50*1,'s');
        self.axM=self.drawAx(50*2,'m');
        self.axL=self.drawAx(50*3,'l');
        self.axR=self.drawAx(50*4,'r');
        
class plotDatum:
    #Class attributes
    w = 5;
    def __init__(self,dwg,x=0,y=0,h=0,stroke="none",fill="black",translate=0):
        self.dwg = dwg
        self.canvas = dwg.canvas
        self.h = h * self.dwg.scalar
        #Coordinates
        self.x = x
        self.y = y
        self.x_svg = self.x+self.w/2
        self.y_svg = self.dwg.y_size - self.y - self.h
        
        #tranformations
        # self.rotAngle = rotate
        # self.x_rotCenter = self.x + self.w/2
        # self.y_rotCenter = self.dwg.y_size - self.y
        # self.rotDot(); #get coordiantes after rotation
        self.translate = translate
        
        #format
        self.stroke = stroke
        self.stroke_width = 4
        self.fill = fill
        
    def addBar(self,):
        #Bar
        self.canvas.add(
            self.canvas.rect(
                insert = (self.x,self.y_svg),
                size = (self.w, self.h),
                stroke_width = self.stroke_width,stroke = self.stroke,fill = self.fill,
                transform='translate({0},0)'.format(self.translate)
            )
        )
    def addDot(self,):
        #Dot
        self.canvas.add(
            self.canvas.circle(
                center = (self.x_svg,self.y_svg),
                r=5,
                stroke_width = self.stroke_width,stroke = self.stroke,fill = self.fill,
                transform='translate({0},0)'.format(self.translate)
            )
        )
        
    def addzDot(self,):
        #Dot for zData (unable to assign to a specific photoreceptor type)
        self.canvas.add(
            self.canvas.circle(
                center = (self.x_svg,self.y_svg),
                r=5,
                stroke_width = self.stroke_width,stroke = self.stroke,fill = self.fill,
                transform='translate({0},0)'.format(self.translate)
            )
        )
        
    def rotDot(self):
        angle = math.radians(self.rotAngle)
        self.x_real = self.x_rotCenter + math.cos(angle) * (self.x_svg - self.x_rotCenter) - math.sin(angle) * (self.y_svg - self.y_rotCenter)
        self.y_real = self.y_rotCenter + math.sin(angle) * (self.x_svg - self.x_rotCenter) + math.cos(angle) * (self.y_svg - self.y_rotCenter)
        

class uBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,"none",dwg.colors.get('uv'),50*0)
        self.x_real = self.x + 50*0
        self.y_real = self.y

class sBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,"none",dwg.colors.get('s'),50*1)
        self.x_real = self.x + 50*1
        self.y_real = self.y

class mBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,"none",dwg.colors.get('m'),50*2)
        self.x_real = self.x + 50*2
        self.y_real = self.y
        
class lBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,"none",dwg.colors.get('l'),50*3)
        self.x_real = self.x + 50*3
        self.y_real = self.y
        
class rBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,"none",dwg.colors.get('r'),50*4)
        self.x_real = self.x + 50*4
        self.y_real = self.y

class zBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center+dwg.z_gap,dwg.y_center,h,"none",dwg.colors.get('z'),0)
        
class plotCell:
    def __init__(self,dwg,cellData):
        self.dwg = dwg
        self.cellData = cellData
        #create structure and gather data
        self.uDat = uBar(self.dwg,self.cellData[0]);
        self.sDat = sBar(self.dwg,self.cellData[1]);
        self.mDat = mBar(self.dwg,self.cellData[2]);
        self.lDat = lBar(self.dwg,self.cellData[3]);
        self.rDat = rBar(self.dwg,self.cellData[4]);
        self.zDat = zBar(self.dwg,self.cellData[5]);
        #plot Polygon
        self.polyCell();
        #plot individual dots 
        self.uDat.addDot();
        self.sDat.addDot();
        self.mDat.addDot();
        self.lDat.addDot();
        self.rDat.addDot();
        self.zDat.addzDot();
        self.uDat.addBar();
        self.sDat.addBar();
        self.mDat.addBar();
        self.lDat.addBar();
        self.rDat.addBar();
        self.zDat.addzDot();
        
    def polyCell(self):
        uStr="M" + str(int(self.uDat.x_real)) + "," + str(int(self.uDat.y_real))
        sStr="L" + str(int(self.sDat.x_real)) + "," + str(int(self.sDat.y_real))
        mStr="L" + str(int(self.mDat.x_real)) + "," + str(int(self.mDat.y_real))
        lStr="L" + str(int(self.lDat.x_real)) + "," + str(int(self.lDat.y_real))
        rStr="L" + str(int(self.rDat.x_real)) + "," + str(int(self.rDat.y_real))
        d_command=" ".join([uStr,sStr,mStr,lStr,rStr, "Z"])
        self.poly = self.dwg.canvas.add(
            svgwrite.path.Path(d=d_command,
            stroke_linecap="round",
            fill="none",stroke = "rgb(0,0,0)",stroke_width=3,stroke_opacity="1")
        )