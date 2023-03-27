using OxyPlot;
using OxyPlot.Axes;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフのX軸設定
    /// </summary>
    internal class XAxisConfig
    {
        /// <summary>
        /// X軸設定を作成する。
        /// </summary>
        /// <param name="min">最小値</param>
        /// <param name="max">最大値</param>
        /// <param name="majorStep">Major Step(MB)</param>
        public XAxisConfig(double min, double max, XAxisMajorStep majorStep)
        {
            Minimum = min;
            Maximum = max;
            MajorStep = majorStep;
        }

        /// <summary>
        /// X軸の最小値を取得する。
        /// </summary>
        internal double Minimum { get; }

        /// <summary>
        /// X軸の最大値を取得する。
        /// </summary>
        internal double Maximum { get; }

        /// <summary>
        /// X軸 Major Step(MB)を取得する。
        /// </summary>
        internal XAxisMajorStep MajorStep { get; }

        /// <summary>
        /// X軸を作成する。
        /// </summary>
        /// <returns>X軸</returns>
        public LinearAxis CreateAxis()
        {
            return new LinearAxis()
            {
                Title = "Position (MB)",
                Position = AxisPosition.Bottom,
                MajorGridlineStyle = LineStyle.Solid,
                MajorGridlineColor = ColorPalette.MajorGridlineColor,
                MajorStep = MajorStep.Value,
                MinorStep = 1,
                MinorGridlineStyle = LineStyle.None,
                Minimum = Minimum,
                Maximum = Maximum,
                TitleFontSize = 20,
                FontSize = 15
            };
        }
    }
}
