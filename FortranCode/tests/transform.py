from lxml import etree
import sys

if sys.version_info.major == 3:
    from io import BytesIO
else:
    from StringIO import StringIO as BytesIO


TAGfile = open(sys.argv[1]+"/Testing/TAG", 'r')
dirname = TAGfile.readline().strip()

xmlfile = open(sys.argv[1]+"/Testing/"+dirname+"/Test.xml", 'rb')
xslfile = open(sys.argv[2], 'rb')

xmlcontent = xmlfile.read()
xslcontent = xslfile.read()

xmldoc = etree.parse(BytesIO(xmlcontent))
xslt_root = etree.XML(xslcontent)
transform = etree.XSLT(xslt_root)

result_tree = transform(xmldoc)
print(result_tree)
