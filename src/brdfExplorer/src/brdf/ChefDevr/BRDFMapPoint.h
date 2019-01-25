#ifndef BRDFPOINT_H
#define BRDFPOINT_H

/**
 * @file BRDFMapPoint.h
 */

#include <QGraphicsItem>
#include <string>
namespace ChefDevr
{
    //TODO : template<bool ui=false>
    /**
    * @brief Point on the BRDFMap corresponding to a BRDF
    */
    class BRDFMapPoint : public QGraphicsItem{
    public:
        explicit BRDFMapPoint (const std::string& name, QColor clr, bool ui);
        ~BRDFMapPoint () override;
    private:
        /**
        * @brief The name of the point
        */
        const std::string _name;

        /**
        * @brief Is point added manually by the user
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
} // namespace ChefDevr

#endif // BRDFPOINT_H
