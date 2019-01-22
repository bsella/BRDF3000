#ifndef BRDFPOINT_H
#define BRDFPOINT_H

#include <QGraphicsItem>
#include <string>

/**
 * Point on the BRDFMap corresponding to a BRDF
 */
//TODO : template<bool ui=false>
class BRDFMapPoint : public QGraphicsItem{
public:
	explicit BRDFMapPoint (const std::string& name, QColor clr, bool ui);
	~BRDFMapPoint () override;
private:
	/**
	 * The name of the point
	 */
	const std::string _name;

	/**
	 * Is point added manually by the user
	 */
	bool userInput;
	/**
	 * @brief Override of base class QGraphicsItem's
	 * fuction : boundingRect
	 * @details Returns a rectangle that defines the
	 * item's boundaries (x; y; width; height)
	 */
	QRectF boundingRect() const override;
	
	/**
	 * @brief Override of base class QGraphicsItem's
	 * fuction : paint
	 * @details This function defines how the item is
	 * drawn on the scene
	 */
	void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;
	
	/**
	 * @brief Triggered when the cursor enters the
	 * bounding rectangle of the item
	 * @details Shows the point's name on the scene
	 * when the cursor is on the point
	 */
	void hoverEnterEvent(QGraphicsSceneHoverEvent*);
	
	/**
	 * @brief Triggered when the cursor leaves the
	 * bounding rectangle of the item
	 * @details Hides the point's name from the scene
	 * when the cursor leaves the point
	 */
	void hoverLeaveEvent(QGraphicsSceneHoverEvent*);
};

#endif // BRDFPOINT_H
