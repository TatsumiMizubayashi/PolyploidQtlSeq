using OxyPlot;
using OxyPlot.Axes;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフのY軸設定
    /// </summary>

    internal class YAxisConfig
    {
        /// <summary>
        /// Y軸設定を作成する。
        /// </summary>
        /// <param name="min">最小値</param>
        /// <param name="max">最大値</param>
        /// <param name="step">Major step</param>
        /// <param name="distance">Axis Title Distance</param>
        public YAxisConfig(double min, double max, double step, double distance = 8)
        {
            Minimum = min;
            Maximum = max;
            MajorStep = step;
            AxisTitleDistance = distance;
        }

        /// <summary>
        /// 最小値を取得する。
        /// </summary>
        internal double Minimum { get; }

        /// <summary>
        /// 最大値を取得する。
        /// </summary>
        internal double Maximum { get; }

        /// <summary>
        /// MajorStepを取得する。
        /// </summary>
        public double MajorStep { get; }

        /// <summary>
        /// タイトルと目盛り間の距離を取得する。
        /// </summary>
        public double AxisTitleDistance { get; }

        /// <summary>
        /// Y軸を作成する。
        /// </summary>
        /// <param name="title">Y軸タイトル</param>
        /// <returns>Y軸</returns>
        public LinearAxis CreateAxis(string title)
        {
            return new LinearAxis()
            {
                Title = title,
                Position = AxisPosition.Left,
                MajorGridlineStyle = LineStyle.Solid,
                MajorGridlineColor = ColorPalette.MajorGridlineColor,
                MajorStep = MajorStep,
                AxisTitleDistance = AxisTitleDistance,
                Minimum = Minimum,
                Maximum = Maximum,
                TitleFontSize = 20,
                FontSize = 15
            };
        }

        /// <summary>
        /// AxisTitleDistanceを補正したYAxisConfigを作成する。
        /// </summary>
        /// <param name="distanceOffset">補正値</param>
        /// <returns>YAxisConfig</returns>
        public YAxisConfig Correct(double distanceOffset)
        {
            return new YAxisConfig(
                Minimum,
                Maximum,
                MajorStep,
                AxisTitleDistance + distanceOffset);
        }
    }
}
