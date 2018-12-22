#ifndef BRDFMAP_H
#define BRDFMAP_H

#include <QOpenGLWidget>
#include <string>

class BRDFMap : public QOpenGLWidget{
public:
	explicit BRDFMap (const std::string& filePath, QWidget* parent = nullptr);
	~BRDFMap () override;
private:
	void paintEvent (QPaintEvent*) override;
};

#endif // BRDFMAP_H
