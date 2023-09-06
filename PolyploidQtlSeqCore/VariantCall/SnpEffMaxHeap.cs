namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEff MaxHeap
    /// </summary>
    internal class SnpEffMaxHeap
    {
        /// <summary>
        /// max heapの最小値
        /// </summary>
        private const int MINIMUM = 2;

        /// <summary>
        /// max heapの最大値
        /// </summary>
        private const int MAXIMUM = 100;

        /// <summary>
        /// SnpEff MaxHeapを作成する。
        /// </summary>
        /// <param name="maxHeap">最大ヒープサイズ(GB)</param>
        public SnpEffMaxHeap(int maxHeap)
        {
            if (maxHeap < MINIMUM || maxHeap > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(maxHeap));

            Value = maxHeap;
        }

        /// <summary>
        /// 最大ヒープサイズ(GB)を取得する。
        /// </summary>
        internal int Value { get; }

    }
}
