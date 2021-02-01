function PaintVector(Point1,Point2,Color)

if size(Point1,1)>size(Point1,2)
    Point1=Point1';
end

if size(Point2,1)>size(Point2,2)
    Point2=Point2';
end


line = [Point2;Point1];

plot3(line(:,1),line(:,2),line(:,3),Color)


end
