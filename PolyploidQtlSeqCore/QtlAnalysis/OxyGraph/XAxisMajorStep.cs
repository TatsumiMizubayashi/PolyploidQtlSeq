namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// X軸の目盛り間隔(MB)
    /// </summary>
    internal class XAxisMajorStep
    {
        /// <summary>
        /// 目盛り間隔の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 目盛り間隔の最大値
        /// </summary>
        private const int MAXIMUM = 50;


        /// <summary>
        /// X軸の目盛り間隔(MB)を作成する。
        /// </summary>
        /// <param name="value">ステップ数(MB)</param>
        public XAxisMajorStep(int value)
        {
            if (value < MINIMUM || value > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(value));

            Value = value;
        }

        /// <summary>
        /// X軸の目盛り間隔(MB)を取得する。
        /// </summary>
        internal int Value { get; }
    }
}
