namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// (全サンプルの）最低Depthしきい値
    /// </summary>
    internal class MinimumDepthThreshold
    {
        /// <summary>
        /// 最低Depthしきい値の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 最低Depthしきい値の最大値
        /// </summary>
        private const int MAXIMUM = 10000;


        /// <summary>
        /// 最低Depthしきい値を作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public MinimumDepthThreshold(int threshold)
        {
            if (threshold < MINIMUM || threshold > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(threshold));

            Value = threshold;
        }


        /// <summary>
        /// 最低Depthしきい値を取得する。
        /// </summary>
        internal int Value { get; }

    }
}
