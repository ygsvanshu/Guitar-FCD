
*************************************************************************************************************************************************************************************************************

GuitarData.h5 - How to use:
The heirarchial structure is as follows

GuitarData.h5
    -> README                   (Dataset: String)
    -> DATA                     (Group)
        -> Row                  (Group) 
            -> Column           (Group)
                -> Day          (Dataset: Integer)
                -> LightPoint   (Dataset: Float [2,1])  
                -> LightLine    (Dataset: Float [2,2])
                -> CameraLine   (Dataset: Float [2,2])
                -> Phi          (Dataset: Float, units: rad)
                -> D            (Dataset: Float, units: mm)
                -> Theta        (Dataset: Float, units: rad)

README:     This information mentioned here
DATA:       The complete data extracted from images
Row:        The row number of the squares drawn on the guitar body
Column:     The column number of the squares drawn on the guitar body
Day:        The day on which this data was recorded; useful to check against image files which are named as per the days
LightPoint: The center point of LED light panel as marked on the calibrated image. Point in format [x0,y0] in mm from the origin (read below for more info)
LightLine:  The line segment parallel to the of LED surface panel as marked on the calibrated image. Points in format [[x1,y1],[x2,y2]] in mm from the origin (read below for more info)
CameraLine: The line segment (usually) parallel (perpendicular in Day5) to the of lens of the camera surface panel as marked on the calibrated image. Points in format [[x1,y1],[x2,y2]] in mm from the origin (read below for more info)
Phi:        The angle subtended by the LightLine with the surface of the guitar (y-axis)
D:          The normal distance between the LightPoint and the guitar surface (y-axis)
Theta:      The angle subtended by the camera view line with the guitar surface normal (x-axis)

The origin is chosen in this manner:
The guitar surface is on the y-axis of the image (x=0)
The LED is always placed at a lower y-coordinate than the camera for any individual data point
After this alignment, the lower left corner of the sheet image is origin (0 mm, 0 mm)

*************************************************************************************************************************************************************************************************************