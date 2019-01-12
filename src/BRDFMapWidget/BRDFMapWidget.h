#ifndef BRDFMAP_H
#define BRDFMAP_H

#include <QGraphicsScene>
#include <string>

class BRDFMapWidget : public QGraphicsScene{
public:
	explicit BRDFMapWidget (const std::string& filePath, QWidget* parent = nullptr);
	~BRDFMapWidget () override;
private:
	void paintEvent (QPaintEvent*) override;
};

#endif // BRDFMAP_H
