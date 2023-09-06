namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// 分布作成時の試行回数
    /// </summary>
    internal class ReplicatesNumber
    {
        /// <summary>
        /// 試行回数の最小値
        /// </summary>
        private const int MINIMUM = 1000;

        /// <summary>
        /// 試行回数の最大値
        /// </summary>
        private const int MAXIMUM = 1_000_000;

        /// <summary>
        /// 分布作成時の試行回数を作成する。
        /// </summary>
        /// <param name="number">試行回数</param>
        public ReplicatesNumber(int number)
        {
            if (number < MINIMUM || number > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(number));

            Value = number;
        }


        /// <summary>
        /// Bulk2個体数を取得する。
        /// </summary>
        internal int Value { get; }

    }
}
