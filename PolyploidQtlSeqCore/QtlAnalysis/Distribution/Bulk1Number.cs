namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// Bulk1個体数
    /// </summary>
    internal class Bulk1Number
    {
        /// <summary>
        /// 個体数の最小値
        /// </summary>
        private const int MINIMUM = 2;

        /// <summary>
        /// 個体数の最大値
        /// </summary>
        private const int MAXIMUM = 1000;


        /// <summary>
        /// Bulk1個体数を作成する。
        /// </summary>
        /// <param name="number">個体数</param>
        public Bulk1Number(int number)
        {
            if (number < MINIMUM || number > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(number));

            Value = number;
        }

        /// <summary>
        /// Bulk1個体数を取得する。
        /// </summary>
        internal int Value { get; }
    }
}
