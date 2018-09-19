import svgwrite
import math

def rotate(x,y,angle,ox,oy):
    angle = math.radians(angle) #convert angle in deg to radians
    x_rot = ox + math.cos(angle) * (x - ox) - math.sin(angle) * (y - oy)
    y_rot = oy + math.sin(angle) * (x - ox) + math.cos(angle) * (y - oy)
    return x_rot,y_rot
    
class hcCanvas:
    filename = 'test-svgwrite'
    x_size = 800
    y_size = 800
    axColor = "rgb(200,200,200)"
    axColorMin = "rgb(230,230,230)"
    axSW=2
    axSWGon=1
    axMax=25
    scalar = axMax*1/3; #to transform from pixels to units
    z_gap = (axMax+15) * scalar
    h_pentagon = (axMax * scalar * (math.sqrt(5+2*math.sqrt(5))/16))
    def __init__(self):
        self.x_center = self.x_size/2
        self.y_center = self.y_size/2
        self.x_origin = self.x_center
        self.y_origin = self.y_size-self.y_center
        self.colors = options = {'uv':'rgb(127,15,126)','s':'rgb(11,36,251)','m':'rgb(15,127,18)','l':'rgb(252,13,27)','r':'rgb(0,0,0)','z':'rgb(253,189,64)'}
        self.canvas = svgwrite.Drawing(filename = self.filename + ".svg",size = (self.x_size, self.y_size))
        self.canvas.add(self.canvas.rect(insert = (0,0),size = (self.x_size, self.y_size),id='frame',stroke="black",fill="rgb(240,240,240)"))
#         self.addOrigin();
        self.addAxis();
        
        
    def addOrigin(self):
        self.origin=plotDatum(self,self.x_center,self.y_center,0,"none",self.axColor);
        self.origin.addDot('origin');
        self.z_origin=plotDatum(self,self.x_center+self.z_gap,self.y_center-(self.axMax * self.scalar),0,"none",self.axColor);
        self.z_origin.addDot('zOrigin');
    
    def drawAx(self,rotation,tag=""): 
        if tag == 'uv': axLabel, axLabelAnchor, axLabelTrans = 'UV', 'middle', 'translate(0,-10)'
        elif tag == 's': axLabel, axLabelAnchor, axLabelTrans = 'S', 'start', 'translate(+5,0)'
        elif tag == 'm': axLabel, axLabelAnchor, axLabelTrans = 'M', 'start', 'translate(2,10)'
        elif tag == 'l': axLabel, axLabelAnchor, axLabelTrans = 'L', 'end', 'translate(-2,10)'
        elif tag == 'r': axLabel, axLabelAnchor, axLabelTrans = 'Rod', 'end', 'translate(-5,0)'
        else: axLabel, axLabelAnchor, axLabelTrans = '', 'start', 'translate(-5,0)'
        # draw axis    
        Ax = self.canvas.add(self.canvas.line(
            start=(self.x_origin,self.y_origin),
            end=(self.x_origin,self.y_origin-self.scalar*self.axMax),
            id='ax',
            stroke=self.axColor,stroke_width=self.axSW,stroke_dasharray=4,
            transform='rotate({0},{1},{2})'.format(rotation,self.x_center,self.y_size-self.y_center)))
        # place label
        self.canvas.add(self.canvas.text(
            axLabel,
            insert=(rotate(self.x_origin, self.y_origin-self.scalar*self.axMax,rotation,self.x_origin,self.y_origin)),
            font_size=20,font_family='sans-serif',text_anchor=axLabelAnchor,
            fill=self.colors.get(tag), alignment_baseline="middle",
            transform=axLabelTrans
        ))
        return Ax
    
    def drawZAx(self):
        zAx = self.canvas.add(self.canvas.line(
            start=(self.x_center+self.z_gap+2.5,self.y_center+(self.axMax * self.scalar)-self.h_pentagon),
            end=(self.x_center+self.z_gap+2.5,self.y_center-self.h_pentagon),
            id='ax',
            stroke=self.axColor,stroke_width=self.axSW,stroke_dasharray=4,
            ))
        self.canvas.add(self.canvas.text(
            '?',
            insert=(self.x_center+self.z_gap+2.5,self.y_center-self.h_pentagon),
            font_size=20,font_family='sans-serif',text_anchor='middle',
            fill=self.colors.get('z'),
            transform='translate(0,-10)'
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
                id='axGon',fill="white",stroke = self.axColor,stroke_width=self.axSWGon)
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
                id='axGonMinor',fill=self.axColorMin,stroke = self.axColor,stroke_width=self.axSWGon)
        )
    
    def drawZAxGon(self,scale):
        AxSquare = self.canvas.add(
            self.canvas.rect(
                insert = (self.x_center+self.z_gap+2.5-20,self.y_center+((self.axMax-scale) * self.scalar)-self.h_pentagon),
                size = (40,scale * self.scalar),
                id='zaxGon',fill="white",stroke = self.axColor,stroke_width=self.axSWGon)
        )
    
    def drawZAxGonMinor(self,scale):
        AxSquare = self.canvas.add(
            self.canvas.rect(
                insert = (self.x_center+self.z_gap+2.5-20,self.y_center+((self.axMax-scale) * self.scalar)-self.h_pentagon),
                size = (40,scale * self.scalar),
                id='zaxGonMinor',fill=self.axColorMin,stroke = self.axColor,stroke_width=self.axSWGon)
        )
        
    
    
    def addAxis(self):
        self.drawAxgonMinor(25);
        self.drawAxgon(20);
        self.drawAxgonMinor(15);
        self.drawAxgon(10);
        self.drawAxgonMinor(5);
        self.drawZAxGonMinor(25);
        self.drawZAxGon(20);
        self.drawZAxGonMinor(15);
        self.drawZAxGon(10);
        self.drawZAxGonMinor(5);
        
#         [self.drawAx(360*facet/5) for facet in range(5)]
        self.drawZAx();
        self.axU=self.drawAx(360*0/5,'uv');
        self.axS=self.drawAx(360*1/5,'s');
        self.axM=self.drawAx(360*2/5,'m');
        self.axL=self.drawAx(360*3/5,'l');
        self.axR=self.drawAx(360*4/5,'r');
        
class plotDatum:
    #Class attributes
    w = 5;
    def __init__(self,dwg,x=0,y=0,h=0,stroke="black",fill="none",rotate=0):
        self.dwg = dwg
        self.canvas = dwg.canvas
        self.h = h * self.dwg.scalar
        #Coordinates
        self.x = x
        self.y = y
        self.x_svg = self.x+self.w/2
        self.y_svg = self.dwg.y_size - self.y - self.h
        
        #tranformations
        self.rotAngle = rotate
        self.x_rotCenter = self.x + self.w/2
        self.y_rotCenter = self.dwg.y_size - self.y
        self.rotDot(); #get coordiantes after rotation
        
        #format
        self.stroke = stroke
        self.stroke_width = 4
        self.fill = fill
        
    def addBar(self,svg_id=''):
        #Bar
        self.canvas.add(
            self.canvas.rect(
                insert = (self.x,self.y_svg),
                size = (self.w, self.h),
                id=svg_id,
                stroke_width = self.stroke_width,stroke = self.stroke,fill = self.fill,
                transform='rotate({0},{1},{2})'.format(self.rotAngle,self.x_rotCenter,self.y_rotCenter)
            )
        )
    def addDot(self,svg_id=''):
        #Dot
        self.canvas.add(
            self.canvas.circle(
                center = (self.x_svg,self.y_svg),
                id=svg_id,
                r=5,
                stroke_width = self.stroke_width,stroke = self.stroke,fill = self.fill,
                transform='rotate({0},{1},{2})'.format(self.rotAngle,self.x_rotCenter,self.y_rotCenter)
            )
        )
        
    def addzDot(self,svg_id=''):
        #Dot for zData (unable to assign to a specific photoreceptor type)
        self.canvas.add(
            self.canvas.circle(
                center = (self.x_svg,self.y_svg+(self.dwg.axMax * self.dwg.scalar)-self.dwg.h_pentagon),
                id=svg_id,
                r=5,
                stroke_width = self.stroke_width,stroke = self.stroke,fill = self.fill,
                transform='rotate({0},{1},{2})'.format(self.rotAngle,self.x_rotCenter,self.y_rotCenter)
            )
        )
        
    def rotDot(self):
        angle = math.radians(self.rotAngle)
        self.x_real = self.x_rotCenter + math.cos(angle) * (self.x_svg - self.x_rotCenter) - math.sin(angle) * (self.y_svg - self.y_rotCenter)
        self.y_real = self.y_rotCenter + math.sin(angle) * (self.x_svg - self.x_rotCenter) + math.cos(angle) * (self.y_svg - self.y_rotCenter)
        

class uBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,dwg.colors.get('uv'),"none",360*0/5)
        # self.addDot('uDat')

class sBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,dwg.colors.get('s'),"none",360*1/5)
        # self.addDot('sDat')

class mBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,dwg.colors.get('m'),"none",360*2/5)
        # self.addDot('mDat')

class lBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,dwg.colors.get('l'),"none",360*3/5)
        # self.addDot('lDat')

class rBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center,dwg.y_center,h,dwg.colors.get('r'),"none",360*4/5)
        # self.addDot('rDat')

class zBar(plotDatum):
    def __init__(self,dwg,h):
        super().__init__(dwg,dwg.x_center+dwg.z_gap,dwg.y_center,h,dwg.colors.get('z'),"none",0)
        self.addzDot('zDat')
        
class plotCell:
    def __init__(self,dwg,cellData):
        self.dwg = dwg
        self.cellData = cellData
        self.uDat = uBar(self.dwg,self.cellData[0]);
        self.sDat = sBar(self.dwg,self.cellData[1]);
        self.mDat = mBar(self.dwg,self.cellData[2]);
        self.lDat = lBar(self.dwg,self.cellData[3]);
        self.rDat = rBar(self.dwg,self.cellData[4]);
        self.zDat = zBar(self.dwg,self.cellData[5]);
        self.polyCell();
        
    def polyCell(self):
        self.poly = self.dwg.canvas.add(
            svgwrite.shapes.Polygon(
                points=[
                    [self.uDat.x_real,self.uDat.y_real],
                    [self.sDat.x_real,self.sDat.y_real],
                    [self.mDat.x_real,self.mDat.y_real],
                    [self.lDat.x_real,self.lDat.y_real],
                    [self.rDat.x_real,self.rDat.y_real]
                ],
                id='cellGon',fill="none",stroke = "rgb(0,0,0)",stroke_width=3,stroke_opacity="0.5")
        )